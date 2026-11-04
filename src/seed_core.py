"""Counted U/E moment evaluation with shared kernels and 80-digit interval arithmetic."""
from fractions import Fraction as F
from decimal import Decimal
from math import comb,log,log1p,exp,sqrt,isfinite
import os,resource,time,json,hashlib
from runtime import check_resources
from pathlib import Path
from interval_core import I,ZERO,ONE,exp_fraction,divide,packed,rdp_width
Z=(0,0,0);NEG=float('-inf')
def validate(s):
    if s.get('bandwidth',s['b'])>s['b'] or s.get('amplitude_dependent') or s.get('heterogeneous') or (s.get('data_mechanism') and s.get('C')=='zero'):raise ValueError('UNSUPPORTED')
    assert all(type(s[k]) is int and s[k]>0 for k in ('T','b','k','a')) and s['b']>=2 and s['a']>=2
    assert 0<F(s['p'])<1 and F(s['sigma2'])>0

def ladd(x,y):
    if x==NEG:return y
    if y==NEG:return x
    hi,lo=max(x,y),min(x,y);return hi+log1p(exp(lo-hi))

class Algebra:
    def __init__(self,cert):self.cert=cert;self.cache={};self.zero=ZERO if cert else NEG;self.one=ONE if cert else 0.
    def add(self,x,y):return x+y if self.cert else ladd(x,y)
    def mul(self,x,y):return x*y if self.cert else x+y
    def prob(self,x):return I.exact(x) if self.cert else log(float(x))
    def ex(self,n):
        if self.cert:
            if type(n) is not int or not 0<=n<=15360:raise ArithmeticError('ENERGY_GRID')
            if n not in self.cache:
                if len(self.cache)>=15361:raise ArithmeticError('EXP_CACHE_CAP')
                self.cache[n]=exp_fraction(F(n,168))
            return self.cache[n]
        return n
    def value(self,x):return packed(x) if self.cert else x

class Work:
    def __init__(self,limit):self.n=0;self.periodic=0;self.boundaries=0;self.cpu=0.;self.limit=limit
    def check(self):
        t=time.process_time();check_resources()
        self.cpu+=time.process_time()-t
    def tick(self):
        self.n+=1
        if self.n>self.limit:raise ArithmeticError('LOGICAL_WORK_CAP')
        if self.n%4096==0:self.periodic+=1;self.check()
    def row(self):self.boundaries+=1;self.check()

def calc(s,arm,cert=False,trace=None):
    validate(s);assert arm in ('U','E')
    a,b,k,T=(s[x] for x in ('a','b','k','T'));p=F(s['p']);sig=F(s['sigma2']);work=Work(2500000 if cert else 25000000);work.check()
    if cert:assert (b,k,p,sig)==(3,2,F(1,64),F(4)) and a<=16 and T<=24
    alg=Algebra(cert);add,mul=alg.add,alg.mul;start=time.process_time()
    q=[F(comb(k,j))*p**j*(1-p)**(k-j) for j in range(k+1)];q=[alg.prob(x) for x in q]
    K=[{Z:alg.one}]
    for h in range(a):
        nxt={}
        for (c,u,v),mass in K[-1].items():
            for j in range(k+1):
                key=(c+(j>0),u+j,v+j*j);work.tick();nxt[key]=add(nxt.get(key,alg.zero),mul(mass,q[j]))
        K.append(nxt)
    pre=time.process_time()-start;D=[2**(b-1-r) for r in range(b)];S=sum(x*x for x in D)
    weights={(Z,)*(b-1):alg.one};rows=[];dp=time.process_time();peak=1
    for t in range(T):
        co=[x if t>=r else 0 for r,x in enumerate(D)]
        if s.get('C',s.get('coefficient_family'))=='zero':co=[0]*b
        if not cert:
            eta=.75 if s.get('C',s.get('coefficient_family'))=='alternating_rows' and (t+1)%2==0 else 1.
            co=[eta*x/sqrt(S) for x in co];den=2*float(sig)
        else:den=168
        counts={'U':0,'A':0,'B':0,'target':0};nxt={}
        def factor(n):return alg.ex(n if cert else n/den)
        if arm=='U':
            for ages,mass in weights.items():
                h=a-sum(z[0] for z in ages)
                for z,prob in K[h].items():
                    groups=(z,)+ages;ss=sum(c*x[1] for c,x in zip(co,groups));vv=sum(c*c*x[2] for c,x in zip(co,groups));en=ss*ss-vv
                    key=(z,)+ages[:-1];val=mul(mul(mass,prob),factor(en));nxt[key]=add(nxt.get(key,alg.zero),val);counts['U']+=1;work.tick()
        else:
            by_y={};g,d=co[0],co[-1]
            for ages,mass in weights.items():
                y=ages[:-1];c,u,v=ages[-1];h=a-sum(z[0] for z in y)-c;sy=sum(co[r+1]*z[1] for r,z in enumerate(y))
                tab=by_y.setdefault(y,{});tab[h,u]=add(tab.get((h,u),alg.zero),mul(mass,factor(d*d*(u*u-v)+2*d*u*sy)));counts['A']+=1;work.tick()
            for y,A in by_y.items():
                sy=sum(co[r+1]*z[1] for r,z in enumerate(y));ey=sy*sy-sum(co[r+1]**2*z[2] for r,z in enumerate(y));B={}
                for (h,u),mass in A.items():
                    for x in range(k*h+1):B[h,x]=add(B.get((h,x),alg.zero),mul(mass,factor(2*g*d*u*x)));counts['B']+=1;work.tick()
                for h in sorted({h for h,u in A}):
                    for z,prob in K[h].items():
                        cp,x,vp=z;key=(z,)+y;val=mul(mul(prob,B[h,x]),factor(ey+2*g*x*sy+g*g*(x*x-vp)));nxt[key]=add(nxt.get(key,alg.zero),val);counts['target']+=1;work.tick()
        weights=nxt;peak=max(peak,len(weights));work.row()
        if k==2 and t>=b-1:
            expected=comb(a+2*b,2*b) if arm=='U' else 3*comb(a+2*b-1,2*b-1)
            assert sum(counts.values())==expected,'FULL_ROW_COUNT'
        row={'row':t+1,'states':len(weights),'counts':counts}
        if trace:
            target=Path(trace)/('row%03d.json'%(t+1));data=[[key,alg.value(value)] for key,value in weights.items()]
            target.write_text(json.dumps(data,allow_nan=False,separators=(',',':')));row.update(trace=str(target),sha256=hashlib.sha256(target.read_bytes()).hexdigest())
        rows.append(row)
    total=alg.zero
    for mass in weights.values():total=add(total,mass)
    if cert:total=I(max(Decimal(1),total.lo),total.hi)
    else:assert isfinite(total) and total>=-1e-10
    return {'status':'PASS','numeric':'outward_interval' if cert else 'ordinary_float_not_certificate','value':alg.value(total),'rows':rows,'peak_states':peak,'preprocess_cpu':pre,'recurrence_cpu':time.process_time()-dp,'monitor':vars(work),'exp_cache_entries':len(alg.cache),'kernel_sizes':[len(x) for x in K]}

def certificate(s,arm):
    if (s['T'],s['b'],s['k'],s['a'],F(s['p']),F(s['sigma2'])) != (24,3,2,16,F(1,64),F(4)) or s.get('C',s.get('coefficient_family','geometric')) != 'geometric':
        raise ValueError('certificate requires the fixed W3 configuration')
    r=calc(s,arm,True);v=I(Decimal(r['value']['lo']),Decimal(r['value']['hi']))
    forward=divide(I.exact(F(15**15,16**16))*v,exp_fraction(15*F(657,1000)))
    reverse=divide(ONE,exp_fraction(F('11.53390625')));width=rdp_width(v,16);target=Decimal('.00001')
    r.update(forward_delta=packed(forward),reverse_delta=packed(reverse),rdp_width=str(width),certificate_status='PASS' if max(forward.hi,reverse.hi)<=target and width<=Decimal('1e-8') else 'UNDECIDED',epsilon='657/1000',delta='1/100000')
    return r
