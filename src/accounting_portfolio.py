"""Deterministic B/U/E privacy-query portfolio and moment-chord conversion."""
import math,time
from workload_core import calc,coefficients,ladd,NEG

def exact_zero(p,zero_C=False):
    assert p==0 or zero_C
    return 0.

def convert(r,a,d):
    return max(0.,r+(-math.log(d)+(a-1)*math.log(a-1)-a*math.log(a))/(a-1))

def best(items):
    assert items and all(math.isfinite(x['epsilon']) and x['epsilon']>=0 for x in items)
    return min(items,key=lambda x:x['epsilon'])

def candidate(label,e,orders=()):return dict(label=label,epsilon=e,orders=list(orders))

def jhs(e,a0,tau):
    if tau==0:return -math.expm1(e-a0) if e<a0 else 0.
    z=(a0-e)/tau
    # Equivalent normal-CDF expression, with stable scaled erfc in the tail.
    t=(tau-z)/math.sqrt(2); A=e-a0+tau*tau/2
    if t>8:
        # No asymptotic approximation: erfc is still representable for our fixed range.
        v=math.erfc(t)
        if v==0: raise ArithmeticError('JENSEN_ERFC_UNDERFLOW')
        second=math.exp(A+math.log(v)-math.log(2))
    else: second=math.exp(A)*.5*math.erfc(t)
    value=.5*math.erfc(-z/math.sqrt(2))-second
    if value < -1e-12 or not math.isfinite(value):raise ArithmeticError('JENSEN_CANCELLATION')
    return max(0.,value)

def root(fun,d):
    if fun(0)<=d:return 0.
    hi=1.
    while fun(hi)>d and hi<1024:hi*=2
    if fun(hi)>d:raise ArithmeticError('ROOT_NOT_BRACKETED')
    lo=0.
    for _ in range(65):
        mid=(lo+hi)/2
        if fun(mid)<=d:hi=mid
        else:lo=mid
    return hi

def screening(triples):
    # Equality from reverse dominance is exact max(P,R)==R; allow it explicitly.
    # Otherwise use a strict margin for the lower-candidate test.
    return all(P<=R or max(P,R)<max(R,L)-1e-8*(1+abs(P)+abs(R)+abs(L)) for P,R,L in triples)

def interpolated(logs,d):
    points=[(1,0.)]+sorted((int(a),f) for a,f in logs.items()); items=[]
    for (l,fl),(h,fh) in zip(points,points[1:]):
        u=(fh-fl)/(h-l);v=fl-u*l;z=u+v-math.log(d)
        star=l if z<=math.log(l) else h if z>=math.log(h) else math.exp(z)
        for a in sorted(set([l,h,star])):
            if a<=1:continue
            # Endpoints use the original log moment; zero-weight neighbours
            # must never receive contribution credit through rounding/tie order.
            if a==l:f,weights=fl,{str(l):1.,str(h):0.}
            elif a==h:f,weights=fh,{str(l):0.,str(h):1.}
            else:f,weights=u*a+v,{str(l):(h-a)/(h-l),str(h):(a-l)/(h-l)}
            item=candidate('moment_chord',convert(f/(a-1),a,d),(l,h))
            item.update(evaluation_order=a,bracket=[l,h],endpoint_weights=weights,
                        active_orders=[x for x in (l,h) if weights[str(x)]>0])
            items.append(item)
    return items

def cheap(spec,logf,deltas):
    T,b,k,p,sig=(spec[x] for x in ('T','b','k','p','sigma'));orders=spec['orders'];co=coefficients(b)
    if not 0<p<1:raise ValueError('general portfolio requires 0<p<1; deterministic QA separate')
    ready=1.;cool=[0.]*(b-1);ex=[];ex2=[]
    for t in range(T):
        ex.append(ready*k*p);ex2.append(ready*(k*p*(1-p)+(k*p)**2))
        pos=ready*(1-(1-p)**k);ready,cool=ready*(1-p)**k+cool[-1],[pos]+cool[:-1]
    norms=[sum(co[r]**2 for r in range(min(b,T-s))) for s in range(T)]
    v=sum(x*y for x,y in zip(ex2,norms));means=[sum(ex[t-r]*co[r] for r in range(min(b,t+1))) for t in range(T)]
    m2=sum(x*x for x in means);a0=v/(2*sig*sig);tau=math.sqrt(m2)/sig;lw=k*T*math.log1p(-p)
    if logf<0:raise ArithmeticError('NEGATIVE_LOGF2')
    c=math.expm1(logf);path={}
    for a in orders:
        z=[0.]*(T+b+1)
        for t in range(T-1,-1,-1):
            terms=[k*math.log1p(-p)+z[t+1]]
            for j in range(1,k+1):
                q=math.comb(k,j)*p**j*(1-p)**(k-j)
                terms.append(math.log(q)+a*(a-1)*j*j*norms[t]/(2*sig*sig)+z[t+b])
            val=NEG
            for term in terms:val=ladd(val,term)
            z[t]=val
        path[a]=z[0]/(a-1)
    out={}
    for d in deltas:
        cs=lambda x:math.log(max(1.,1+x/(4*d)-d))
        pos=[candidate('F2_CS',cs(c),(2,)),candidate('F2_RDP',convert(logf,2,d),(2,))]
        pos.extend(candidate('path',convert(path[a],a,d),(a,)) for a in orders)
        rev=[candidate('reverse_CS',cs(c*math.exp(-lw)),(2,)),candidate('zero_HS',max(0.,math.log1p(-d)-lw)),candidate('Jensen_HS',root(lambda e:jhs(e,a0,tau),d))]
        for a in orders:
            lk=NEG
            for j in range(1,a):lk=ladd(lk,math.log(a-j)-j*lw)
            lc=math.log(c)+lk if c else NEG
            chi=ladd(0.,lc)/(a-1)
            r=min(-lw,(v+(a-1)*m2)/(2*sig*sig),chi)
            rev.append(candidate('reverse_RDP',convert(r,a,d),(a,)))
        P,R=best(pos),best(rev)
        out[str(d)]={'forward':P,'reverse':R,'epsilon':max(P['epsilon'],R['epsilon']),'lower_missing':convert(logf,max(orders),d)}
    return out,dict(v=v,m_norm2=m2,log_w0=lw,paths=path,single_chain_steps=T,path_steps=T*len(orders),root_iterations=65*len(deltas))

def workload(spec,arm,deltas):
    assert arm in ('B','U','E');cache={};start=time.process_time()
    first=calc(dict(spec,a=2,arm='U'),cache);logs={2:first['logF']};details={'2':first}
    cb=time.process_time();base,common=cheap(spec,logs[2],deltas);cheap_cpu=time.process_time()-cb
    skip=screening([(x['forward']['epsilon'],x['reverse']['epsilon'],x['lower_missing']) for x in base.values()])
    # The common screen cannot hide an F1/F2 chord gain absent from B.
    skip=skip and all(max(min([x['forward']['epsilon']]+[c['epsilon'] for c in interpolated(logs,float(d))]),x['reverse']['epsilon'])==x['epsilon'] for d,x in base.items())
    high_start=time.process_time()
    if arm!='B' and not skip:
        for a in spec['orders'][1:]:
            d=calc(dict(spec,a=a,arm=arm),cache);logs[a]=d['logF'];details[str(a)]=d
    high_cpu=time.process_time()-high_start;conv=time.process_time();out={}
    for delta in deltas:
        s=base[str(delta)];options=[s['forward']]
        if arm!='B':options+=interpolated(logs,delta)
        P=best(options);R=s['reverse'];eps=max(P['epsilon'],R['epsilon'])
        if arm!='B' and skip:
            assert abs(eps-s['epsilon'])<=1e-8*(1+abs(s['epsilon']))
        out[str(delta)]={'forward':P,'reverse':R,'epsilon':eps,'baseline_epsilon':s['epsilon']}
    for l,h in zip(deltas,deltas[1:]):assert out[str(h)]['epsilon']+1e-8>=out[str(l)]['epsilon']
    return {'schema':'ic-workload-result-v1','numeric_status':'ordinary_float_theoretical_bounds_not_certificate','arm':arm,'spec':spec,'screen':skip,'moment_calls':len(logs),'moments':details,'queries':out,'base':base,'common':common,'kernel_levels_built':len(cache['logs']),'cheap_cpu':cheap_cpu,'high_cpu':high_cpu,'conversion_cpu':time.process_time()-conv,'process_cpu_before_serialization':time.process_time()-start,'status':'PASS'}
