"""Primitive, independent-oracle, support-count and input-validation checks."""
from fractions import Fraction as F
from decimal import Decimal
import math
from seed_core import calc,validate
import path_oracle

def near(x,y):return abs(x-y)<=1e-10+1e-9*abs(y)
def qa(q,trace_root):
    if q['id']=='A1':
        from certificate_qa import primitive
        result=primitive()
        from interval_core import exp_fraction
        for invalid in ('1',.1):
            try:exp_fraction(invalid)
            except (TypeError,ValueError,AttributeError):pass
            else:raise AssertionError('NON_RATIONAL_EXP_ACCEPTED')
        return result
    if q['id']=='A2':
        from exact_encoding import encode,decode,roundtrip_qa
        encoding_qa=roundtrip_qa()
        s=q;lo,hi,n=path_oracle.rational(s);out={}
        assert decode(encode(lo))==lo and decode(encode(hi))==hi
        for arm in ('U','E'):
            target=trace_root/arm;target.mkdir();r=calc(s,arm,True,target);v=r['value'];l,h=F(Decimal(v['lo'])),F(Decimal(v['hi']))
            if h<lo or l>hi:raise ArithmeticError('ORACLE_DISJOINT')
            if not l<=lo<=hi<=h:return {'status':'UNDECIDED','reason':'ORACLE_CONTAINMENT_UNDECIDED','result':r}
            out[arm]=r
        return {'status':'PASS','results':out,'oracle':{'lo':encode(lo),'hi':encode(hi),'terms':n,'roundtrip_qa':encoding_qa}}
    if q['id']=='A3':
        out=[]
        for s in q['cases']:
            references={str(p):path_oracle.numeric(s,p) for p in (80,120)};assert near(float(references['80'][0]),float(references['120'][0]));r={}
            for arm in ('U','E'):
                dest=trace_root/(s['id']+arm);dest.mkdir();v=calc(s,arm,False,dest);assert near(v['value'],float(references['120'][0]));r[arm]=v
            out.append({'spec':s,'references':references,'results':r})
        return {'status':'PASS','cases':out}
    assert q['id']=='A4';out={}
    for arm in ('U','E'):
        dest=trace_root/arm;dest.mkdir();out[arm]=calc(q,arm,False,dest);assert near(out[arm]['value'],0.)
    rejected=[]
    for extra in ({'bandwidth':5},{'amplitude_dependent':True},{'heterogeneous':True},{'data_mechanism':True}):
        try:validate(dict(q,**extra))
        except ValueError as e:assert str(e)=='UNSUPPORTED';rejected.append(extra)
        else:raise AssertionError('INVALID_MODEL_ACCEPTED')
    return {'status':'PASS','results':out,'rejected':rejected}
