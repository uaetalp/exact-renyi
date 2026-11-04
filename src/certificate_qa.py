"""Independent rational Taylor fixtures for the interval primitives."""
from fractions import Fraction as F
from decimal import Decimal
from itertools import product
from math import factorial
from interval_core import I,ONE,exp_fraction,divide,rdp_width,packed

def taylor(x):
    assert F(0)<=x<=1
    n=200;s=sum((x**j/F(factorial(j)) for j in range(n+1)),F(0))
    tail=x**(n+1)/F(factorial(n+1))/(1-x/F(n+2))
    return s,s+tail

def independent_exp(x):
    if x<=1:return taylor(x)
    assert x<=16
    lo,hi=taylor(x/16)
    for _ in range(4):lo,hi=lo*lo,hi*hi
    return lo,hi

def contains(interval,lo,hi):assert F(interval.lo)<=lo<=hi<=F(interval.hi)

def primitive():
    contains(I.exact(F(1,3))+I.exact(F(1,7)),F(10,21),F(10,21))
    contains(I.exact(F(1,3))*I.exact(F(1,7)),F(1,21),F(1,21))
    contains(divide(I.exact(F(1,3)),I.exact(F(1,7))),F(7,3),F(7,3))
    out=ONE
    for _ in range(100):out=out+I.exact(F(1,10**81))
    contains(out,1+F(100,10**81),1+F(100,10**81))
    for x in [F(0),F(1,3),F(2),15*F(657,1000),F('11.53390625')]:contains(exp_fraction(x),*independent_exp(x))
    assert exp_fraction(F(0))==ONE and rdp_width(ONE,16)==0
    for fn in [lambda:exp_fraction(F(-1)),lambda:divide(ONE,I.exact(0)),lambda:I(Decimal('NaN'),Decimal(1)),lambda:I(Decimal(2),Decimal(1))]:
        try:fn()
        except (ValueError,ArithmeticError):pass
        else:raise AssertionError('invalid input accepted')
    return dict(status='PASS',moment_calls=0,updates=0,reference='independent rational Taylor N200 plus four squares')

