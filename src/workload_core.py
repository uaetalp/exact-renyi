"""Earlier floating-point U/E implementation with shared kernel caching."""
import hashlib, json, math, os, resource, time
from runtime import check_resources, peak_rss_bytes

NEG = float('-inf')
ZERO = (0, 0, 0)


def ladd(x, y):
    if x == NEG: return y
    if y == NEG: return x
    hi, lo = max(x, y), min(x, y)
    return hi + math.log1p(math.exp(lo - hi))


class Monitor:
    def __init__(self): self.updates = self.periodic = self.boundaries = 0; self.cpu = 0.0
    def check(self):
        t = time.process_time()
        check_resources()
        self.cpu += time.process_time() - t
    def tick(self):
        self.updates += 1
        if self.updates % 4096 == 0: self.periodic += 1; self.check()
    def boundary(self): self.boundaries += 1; self.check()
    def record(self): return vars(self).copy()


def kernels(a, k, p, mon, cache=None):
    q = [math.comb(k, j) * p**j * (1-p)**(k-j) for j in range(k+1)]
    if cache is None: cache = {}
    assert cache.get("identity", (k,p)) == (k,p)
    cache["identity"]=(k,p)
    tables=cache.setdefault("raw", [{ZERO:1.0}])
    logs=cache.setdefault("logs", [{ZERO:0.0}])
    for _ in range(len(tables)-1,a):
        nxt = {}
        for (c, u, v), w in tables[-1].items():
            for j, prob in enumerate(q):
                if prob == 0: continue
                key = (c+(j > 0), u+j, v+j*j)
                nxt[key] = nxt.get(key, 0.0) + w*prob; mon.tick()
        assert all(w > 0 and math.isfinite(w) for w in nxt.values())
        tables.append(nxt); logs.append({key:math.log(w) for key,w in nxt.items()}); mon.boundary()
    return logs[:a+1]


def coefficients(b, ident=False):
    if ident: return (1.0,) + (0.0,)*(b-1)
    den = math.sqrt(sum(4.0**(-j) for j in range(b)))
    return tuple(2.0**(-j)/den for j in range(b))


def energy(groups, coeff, sigma):
    s = sum(x[1]*d for x, d in zip(groups, coeff))
    selfsum = sum(x[2]*d*d for x, d in zip(groups, coeff))
    return (s*s-selfsum)/(2*sigma*sigma)


def calc(spec, cache=None, full=False):
    mon = Monitor(); mon.boundary()
    t0 = time.process_time()
    a, k, b, T, p, sigma = (spec[x] for x in ('a', 'k', 'b', 'T', 'p', 'sigma'))
    coeff = coefficients(b, spec.get('identity_C', False))
    before = len(cache.get("logs", [])) if cache is not None else 0
    kw = kernels(a, k, p, mon, cache)
    pre = time.process_time()-t0; dp = time.process_time()
    initial = tuple(tuple(z) for z in spec.get('initial_ages', [ZERO]*(b-1)))
    weights = {initial: 0.0}; rows = []; peak = 1
    for t in range(T):
        nextw = {}; count = {'U_edges': 0, 'A': 0, 'B': 0, 'target': 0}
        active = tuple(coeff[r] if t >= r else 0.0 for r in range(b))
        if spec['arm'] == 'U':
            for ages, oldw in weights.items():
                h = a-sum(z[0] for z in ages)
                for z, lp in kw[h].items():
                    key = (z,) + ages[:-1]
                    val = oldw+lp+energy((z,)+ages, active, sigma)
                    nextw[key] = ladd(nextw.get(key, NEG), val)
                    count['U_edges'] += 1; mon.tick()
            mon.boundary()
        else:
            # Keep all younger ages. The last old group alone expires here.
            by_y = {}
            for ages, oldw in weights.items():
                y, (c, u, v) = ages[:-1], ages[-1]
                n = a-sum(z[0] for z in y); h = n-c
                d = active[-1]; sy = sum(active[r+1]*z[1] for r, z in enumerate(y))
                val = oldw+(d*d*(u*u-v)+2*d*u*sy)/(2*sigma*sigma)
                tab = by_y.setdefault(y, {})
                tab[h,u] = ladd(tab.get((h,u), NEG), val)
                count['A'] += 1; mon.tick()
            mon.boundary()
            for y, A in by_y.items():
                sy = sum(active[r+1]*z[1] for r,z in enumerate(y))
                ey = energy(y, active[1:-1], sigma)
                d, g = active[-1], active[0]
                B = {}
                for (h,u), val in A.items():
                    for x in range(k*h+1):
                        key = h,x
                        B[key] = ladd(B.get(key,NEG), val+g*d*u*x/(sigma*sigma))
                        count['B'] += 1; mon.tick()
                mon.boundary()
                for h in sorted({h for h,u in A}):
                    for (cp,x,vp), lp in kw[h].items():
                        z = (cp,x,vp); key = (z,)+y
                        val = lp+B[h,x]+ey+g*x*sy/(sigma*sigma)+g*g*(x*x-vp)/(2*sigma*sigma)
                        nextw[key] = ladd(nextw.get(key,NEG), val)
                        count['target'] += 1; mon.tick()
                mon.boundary()
        weights = nextw; peak = max(peak,len(weights)); mon.boundary()
        # Full per-row weights allow tolerant, per-key comparison across serial workers.
        entries = [[list(z) for z in key] + [val] for key,val in sorted(weights.items())] if full else None
        rows.append({'row':t, 'states':len(weights), **count, 'weights':entries})
    logF = NEG
    for val in weights.values(): logF = ladd(logF,val)
    assert math.isfinite(logF)
    return {'logF':logF, 'peak_states':peak,'rows':rows, 'monitor':mon.record(),
            'preprocess_cpu':pre,'recurrence_cpu':time.process_time()-dp,
            'worker_peak_rss_bytes':peak_rss_bytes(),
            'kernel_levels_generated':[int(h>=before) for h in range(a+1)]}


if __name__ == '__main__':
    import sys
    assert len(sys.argv) == 2
    print(json.dumps(calc(json.loads(sys.argv[1])), allow_nan=False, separators=(',',':')))
