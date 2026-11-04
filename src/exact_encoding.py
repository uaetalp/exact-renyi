"""Lossless tagged Fraction hexadecimal encoding."""
from fractions import Fraction
import sys

def encode(x):
    if not isinstance(x,Fraction):raise TypeError('Fraction required')
    return {'encoding':'fraction-hex-v1','numerator_hex':hex(x.numerator),'denominator_hex':hex(x.denominator)}

def decode(x):
    if x.get('encoding')!='fraction-hex-v1':raise ValueError('unknown encoding')
    n=int(x['numerator_hex'],16);d=int(x['denominator_hex'],16)
    if d<=0:raise ValueError('positive denominator required')
    return Fraction(n,d)

def roundtrip_qa():
    before=sys.get_int_max_str_digits() if hasattr(sys,'get_int_max_str_digits') else None
    values=[Fraction(0),Fraction(1,3),Fraction(-7,11),Fraction((1<<20000)+1,(1<<10003)+3)]
    for v in values:
        assert decode(encode(v))==v
        assert decode(encode(-v))==-v
    after=sys.get_int_max_str_digits() if hasattr(sys,'get_int_max_str_digits') else None
    assert before==after
    return {'status':'PASS','encoding':'fraction-hex-v1','fixtures':8,'integer_limit_before':before,'integer_limit_after':after}
