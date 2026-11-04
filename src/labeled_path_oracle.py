"""Independent labelled-path oracle for the earlier small fixtures."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
import json


def paths(T, b, k, p):
    q = [F(__import__('math').comb(k,j))*p**j*(1-p)**(k-j) for j in range(k+1)]
    def visit(t, ready, seq, mass):
        if t == T:
            yield tuple(seq), mass; return
        if ready > t:
            yield from visit(t+1, ready, seq+[0], mass)
        else:
            for j in range(k+1):
                yield from visit(t+1, t+b if j else t+1, seq+[j], mass*q[j])
    return list(visit(0,0,[],F(1)))


def gram_energy(left, right, b, sigma):
    # C[t,t-r]=2^-r/sqrt(Z). The shared sqrt cancels in the inner product.
    Z = sum((F(1,4**j) for j in range(b)),F(0))
    def raw(x,t): return sum((F(x[t-r],2**r) for r in range(b) if t-r>=0),F(0))
    return sum((raw(left,t)*raw(right,t)/Z for t in range(len(left))),F(0))/sigma**2


def value(T,b,k,p,sigma,precision,a=3):
    ps=paths(T,b,k,p)
    with localcontext() as ctx:
        ctx.prec=precision; total=Decimal(0)
        for draws in product(ps,repeat=a):
            prob=F(1); energy=F(0)
            for i,(x,px) in enumerate(draws):
                prob*=px
                for y,py in draws[:i]: energy+=gram_energy(x,y,b,sigma)
            total+=(Decimal(prob.numerator)/Decimal(prob.denominator))*(Decimal(energy.numerator)/Decimal(energy.denominator)).exp()
        return {'paths':len(ps),'labelled_terms':len(ps)**a,'logF':str(total.ln())}

def frozen():
    out={}
    for b,n in [(2,11),(3,7)]:
        x=value(3,b,2,F(1,8),F(2),80); y=value(3,b,2,F(1,8),F(2),120)
        assert x['paths']==y['paths']==n and x['labelled_terms']==n**3
        assert abs(Decimal(x['logF'])-Decimal(y['logF']))<=Decimal('1e-25')
        out[str(b)]={'d80':x,'d120':y}
    return out
