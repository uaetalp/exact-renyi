"""Independent labelled-trajectory oracle using rational and high-precision arithmetic."""
from fractions import Fraction as F
from itertools import product
from math import comb,factorial
from decimal import Decimal,localcontext

def trajectories(s):
    T,b,k=(s[x] for x in ('T','b','k'));p=F(s['p']);out=[]
    def visit(t,ready,x,prob):
        if t==T:out.append((tuple(x),prob));return
        for j in ([0] if t<ready else range(k+1)):
            q=F(1) if t<ready else comb(k,j)*p**j*(1-p)**(k-j)
            visit(t+1,t+b if j else ready,x+[j],prob*q)
    visit(0,0,[],F(1));assert sum(p for x,p in out)==1;return out

def terms(s):
    ps=trajectories(s);a,b,T=(s[x] for x in ('a','b','T'));assert len(ps)**a<=10000
    D=[2**(b-1-r) for r in range(b)];S=sum(x*x for x in D);cache={}
    for x,p in ps:
        cache[x]=tuple(sum((F(x[t-r]*D[r]) for r in range(min(t+1,b))),F(0))*(F(3,4) if s.get('C')=='alternating_rows' and (t+1)%2==0 else 1) for t in range(T))
    energies={}
    for seq in product(ps,repeat=a):
        probability=F(1)
        for x,p in seq:probability*=p
        en=sum((cache[seq[i][0]][t]*cache[seq[j][0]][t] for i in range(a) for j in range(i+1,a) for t in range(T)),F(0))/(S*F(s['sigma2']))
        energies[en]=energies.get(en,F(0))+probability
    return energies,len(ps)**a

def rational_exp(x):
    assert 0<=x<=16
    y=x/16;n=200;term=F(1);lo=term
    for j in range(1,n+1):term=term*y/j;lo+=term
    hi=lo+term*y/(n+1)/(1-y/(n+2))
    for _ in range(4):lo,hi=lo*lo,hi*hi
    return lo,hi

def rational(s):
    es,n=terms(s);lo=F(0);hi=F(0)
    for en,p in es.items():l,h=rational_exp(en);lo+=p*l;hi+=p*h
    return lo,hi,n

def numeric(s,precision):
    es,n=terms(s)
    with localcontext() as c:
        c.prec=precision;out=Decimal(0)
        for en,p in es.items():out+=(Decimal(p.numerator)/Decimal(p.denominator))*(Decimal(en.numerator)/Decimal(en.denominator)).exp()
        return str(out.ln()),n
