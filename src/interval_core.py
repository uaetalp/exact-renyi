"""Outward Decimal arithmetic for nonnegative quantities and exact rational inputs."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import (
    Decimal, Context, InvalidOperation, Overflow, ROUND_CEILING,
    ROUND_FLOOR, ROUND_HALF_EVEN, Subnormal, Underflow, localcontext,
)
from fractions import Fraction
from math import comb, factorial
from typing import Iterator


PRECISION = 80
TRAPS = {InvalidOperation: True, Overflow: True, Underflow: True,
         Subnormal: True}


def context(rounding: str) -> Context:
    c = Context(prec=PRECISION, rounding=rounding,
                Emin=-999999, Emax=999999)
    for signal, enabled in TRAPS.items():
        c.traps[signal] = enabled
    return c


DOWN = context(ROUND_FLOOR)
UP = context(ROUND_CEILING)
NEAR = context(ROUND_HALF_EVEN)


def checked(value: Decimal) -> Decimal:
    if not value.is_finite() or value.is_subnormal():
        raise ArithmeticError("non-normal/non-finite interval endpoint")
    return value


def rational_endpoint(x: Fraction, direction: str) -> Decimal:
    if direction not in ("down", "up"):
        raise ValueError("endpoint direction must be down or up")
    if not isinstance(x, Fraction):
        raise TypeError("endpoint input must be an exact Fraction")
    c = DOWN if direction == "down" else UP
    with localcontext(c):
        return checked(Decimal(x.numerator) / Decimal(x.denominator))


@dataclass(frozen=True)
class I:
    lo: Decimal
    hi: Decimal

    def __post_init__(self) -> None:
        checked(self.lo)
        checked(self.hi)
        if self.lo > self.hi or self.lo < 0:
            raise ArithmeticError("invalid nonnegative interval")

    @classmethod
    def exact(cls, x: int | Fraction) -> I:
        f = Fraction(x)
        if f < 0:
            raise ValueError("negative nonnegative quantity")
        return cls(rational_endpoint(f, "down"),
                   rational_endpoint(f, "up"))

    def __add__(self, other: I) -> I:
        with localcontext(DOWN):
            low = checked(self.lo + other.lo)
        with localcontext(UP):
            high = checked(self.hi + other.hi)
        return I(low, high)

    def __mul__(self, other: I) -> I:
        with localcontext(DOWN):
            low = checked(self.lo * other.lo)
        with localcontext(UP):
            high = checked(self.hi * other.hi)
        return I(low, high)

    def power(self, n: int) -> I:
        if n < 0:
            raise ValueError("negative power")
        out = I.exact(1)
        for _ in range(n):
            out = out * self
        return out


ZERO = I.exact(0)
ONE = I.exact(1)


def exp_rational_interval(low: Fraction, high: Fraction) -> I:
    """Monotone exp on rational bounds, nearest plus adjacent ULPs."""
    if low < 0 or low > high:
        raise ValueError("cross energy must be nonnegative")
    if high == 0:
        return ONE
    arg = I(rational_endpoint(low, "down"), rational_endpoint(high, "up"))
    with localcontext(NEAR):
        a, b = checked(arg.lo.exp()), checked(arg.hi.exp())
    with localcontext(NEAR):
        return I(checked(a.next_minus()), checked(b.next_plus()))


def exp_fraction(x: Fraction) -> I:
    return exp_rational_interval(x, x)


def log2_rational_bounds(terms: int = 120) -> tuple[Fraction, Fraction]:
    """2*atanh(1/3), with an exact positive geometric remainder bound."""
    if terms < 1:
        raise ValueError("need at least one series term")
    partial = 2*sum((Fraction(1, 3**(2*n+1)*(2*n+1))
                     for n in range(terms)), Fraction(0))
    first_omitted = 2*Fraction(1, 3**(2*terms+1)*(2*terms+1))
    return partial, partial + first_omitted/Fraction(8, 9)


def log2_variance_bounds(terms: int = 120) -> tuple[Fraction, Fraction]:
    """The ideal sigma²=1/(2 log 2), contained by rational endpoints."""
    low, high = log2_rational_bounds(terms)
    return Fraction(1, 2)/high, Fraction(1, 2)/low


def product_lower(a: Decimal, b: Decimal) -> Decimal:
    with localcontext(DOWN):
        return checked(a*b)


def product_upper(a: Decimal, b: Decimal) -> Decimal:
    with localcontext(UP):
        return checked(a*b)


def rdp_width(moment: I, a: int) -> Decimal:
    if type(a) is not int or a < 2 or moment.lo < 1:
        raise ArithmeticError("invalid moment for RDP width")
    if moment.lo == moment.hi:
        return Decimal(0)
    with localcontext(NEAR):
        low_log = Decimal(0) if moment.lo == 1 else checked(moment.lo.ln()).next_minus()
        high_log = Decimal(0) if moment.hi == 1 else checked(moment.hi.ln()).next_plus()
        low_log, high_log = checked(low_log), checked(high_log)
    with localcontext(UP):
        return checked((high_log - low_log) / Decimal(a - 1))



def divide(x: I, y: I) -> I:
    if y.lo <= 0:raise ArithmeticError("positive denominator required")
    with localcontext(DOWN):lo=checked(x.lo/y.hi)
    with localcontext(UP):hi=checked(x.hi/y.lo)
    return I(lo,hi)

def packed(x: I):return {"lo":str(x.lo),"hi":str(x.hi)}
