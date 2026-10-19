"""Exact transcription of the uploaded AFH supplement, lines 803–922.

Development/regeneration only: SymPy is NOT a production dependency. Functions
return exp(y)*sqrt(4*pi) times the supplement's reduced single-particle element.
Radial hypergeometrics terminate algebraically; no interpolation/fitting/chop.
All source decimal OBDMEs are interpreted as their written decimal rationals.
"""
from functools import lru_cache

import sympy as s
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j

y = s.Symbol('y', nonnegative=True)
half = s.Rational(1, 2)


def phase(power):
    power = s.sympify(power)
    if not power.is_Integer:
        raise ValueError(f'noninteger phase exponent {power}')
    return s.Integer(-1)**power


def qnorm(j):
    return s.sqrt(2*j+1)


@lru_cache(None)
def three(a, b, c):
    try:
        return wigner_3j(a, b, c, 0, 0, 0)
    except ValueError:
        return s.S.Zero


@lru_cache(None)
def six(a, b, c, d, e, f):
    try:
        return wigner_6j(a, b, c, d, e, f)
    except ValueError:
        return s.S.Zero


@lru_cache(None)
def nine(a, b, c, d, e, f, g, h, i):
    try:
        return wigner_9j(a, b, c, d, e, f, g, h, i)
    except ValueError:
        return s.S.Zero


def orbit(N, two_j):
    j = s.Rational(two_j, 2)
    l = int(j+half) if (N-int(j+half)) % 2 == 0 else int(j-half)
    n = (N-l)//2+1
    if n <= 0 or 2*(n-1)+l != N:
        raise ValueError('unphysical oscillator orbit')
    return n, l, j


def hyp(a, b):
    a, b = s.sympify(a), s.sympify(b)
    if not a.is_Integer or a > 0:
        raise ValueError(f'hypergeometric does not terminate: {a}, {b}')
    return sum((s.rf(a, k)/s.rf(b, k))*y**k/s.factorial(k) for k in range(1-int(a)))


@lru_cache(None)
def radial(np, lp, n, l, L, derivative=0):
    """BesselElement (0), Minus (-1), Plus (+1), with exp(-y) removed."""
    if L < 0:
        return s.S.Zero
    prefactor = (s.Integer(2)**(L-(derivative != 0)) / s.factorial2(2*L+1)
        * y**s.Rational(L-(derivative != 0), 2)
        * s.sqrt(s.factorial(np-1)*s.factorial(n-1)
                 * s.gamma(np+lp+half)*s.gamma(n+l+half)))
    total = 0
    for m in range(n):
        for mp in range(np):
            common = phase(m+mp)/(s.factorial(m)*s.factorial(mp)
                *s.factorial(n-1-m)*s.factorial(np-1-mp)
                *s.gamma(l+m+s.Rational(3, 2))*s.gamma(lp+mp+s.Rational(3, 2)))
            if derivative == 0:
                term = (s.gamma(s.Rational(l+lp+L+2*m+2*mp+3, 2))
                        * hyp(s.Rational(L-lp-l-2*mp-2*m, 2), L+s.Rational(3, 2)))
            else:
                term = -s.Rational(l+lp+L+2*m+2*mp+2, 2)*hyp(
                    s.Rational(L-lp-l-2*mp-2*m-1, 2), L+s.Rational(3, 2))
                multiplier = 2*m if derivative == -1 else 2*l+2*m+1
                if multiplier:
                    term += multiplier*hyp(s.Rational(L-lp-l-2*mp-2*m+1, 2), L+s.Rational(3, 2))
                term *= s.gamma(s.Rational(l+lp+L+2*m+2*mp+2, 2))
            total += common*term
    return s.expand(prefactor*total)


def sigma_l(np, lp, jp, n, l, j, K, L):
    angular = (phase(lp)*s.sqrt((2*j+1)*(2*jp+1)*(2*l+1)*(2*lp+1)*(2*K+1)*(2*L+1)*6)
               *three(lp, L, l)*nine(lp, l, L, half, half, 1, jp, j, K))
    return s.S.Zero if angular == 0 else angular*radial(np, lp, n, l, L)


def phi_piece(np, lp, jp, n, l, j, K, upper, direction):
    L = K+1 if upper else K-1
    if L < 0 or (direction == -1 and l == 0):
        return s.S.Zero
    ell = l+direction
    angular = qnorm(ell)*s.sqrt(l+1 if direction == 1 else l)*three(lp, L, ell)
    if angular == 0:
        return s.S.Zero
    indices = range(K, K+2) if upper else range(K-1, K+1)
    angular *= sum(phase(K-u+(direction == 1))*(2*u+1)
        *six(L, 1, u, 1, K, 1)*six(L, 1, u, l, lp, ell)
        *nine(lp, l, u, half, half, 1, jp, j, K) for u in indices)
    return (s.S.Zero if angular == 0 else angular
            *radial(np, lp, n, l, L, -1 if direction == 1 else 1))


@lru_cache(None)
def element(operator, Np, two_jp, N, two_j, K):
    """PhysicalConditionsNormal/Abnormal and all six elastic operators."""
    np, lp, jp = orbit(Np, two_jp)
    n, l, j = orbit(N, two_j)
    abnormal = operator in ('Sigma_prime', 'Sigma_double_prime', 'Delta')
    if abs(j-jp) > K or K > j+jp or (lp+l+K+abnormal) % 2:
        return s.S.Zero
    if operator == 'M':
        angular = (phase(half+j+K)*s.sqrt((2*j+1)*(2*jp+1)*(2*l+1)*(2*lp+1)*(2*K+1))
                   *three(lp, K, l)*six(lp, jp, half, j, l, K))
        value = s.S.Zero if angular == 0 else angular*radial(np, lp, n, l, K)
    elif operator in ('Phi_double_prime', 'Phi_tilde_prime'):
        upper = sum(phi_piece(np, lp, jp, n, l, j, K, True, d) for d in (1, -1))
        lower = 0 if K == 0 else sum(phi_piece(np, lp, jp, n, l, j, K, False, d) for d in (1, -1))
        overall = phase(lp+1)*6*qnorm(lp)*qnorm(jp)*qnorm(j)
        if operator == 'Phi_double_prime':
            value = overall*(qnorm(K+1)*s.sqrt(K+1)*upper
                              +(qnorm(K-1)*s.sqrt(K)*lower if K else 0))
        else:
            value = overall*(-qnorm(K+1)*s.sqrt(K)*upper
                              +(qnorm(K-1)*s.sqrt(K+1)*lower if K else 0))
            value += sigma_l(np, lp, jp, n, l, j, K, K)/2
    elif operator == 'Delta':
        angular = phase(K+j+half)*qnorm(lp)*qnorm(jp)*qnorm(j)*qnorm(K)**2*six(lp, jp, half, j, l, K)
        total = 0
        for direction in (1, -1):
            if direction == -1 and l == 0:
                continue
            ell = l+direction
            part = ((-1 if direction == 1 else 1)*s.sqrt(l+1 if direction == 1 else l)
                    *qnorm(ell)*six(K, 1, K, l, lp, ell)*three(lp, K, ell))
            if part != 0:
                total += part*radial(np, lp, n, l, K, -1 if direction == 1 else 1)
        value = angular*total
    elif operator == 'Sigma_prime':
        value = (-s.sqrt(K)*sigma_l(np, lp, jp, n, l, j, K, K+1)
                 +(s.sqrt(K+1)*sigma_l(np, lp, jp, n, l, j, K, K-1) if K else 0))/qnorm(K)
    elif operator == 'Sigma_double_prime':
        value = (s.sqrt(K+1)*sigma_l(np, lp, jp, n, l, j, K, K+1)
                 +(s.sqrt(K)*sigma_l(np, lp, jp, n, l, j, K, K-1) if K else 0))/qnorm(K)
    else:
        raise ValueError(operator)
    return s.simplify(s.expand(value))
