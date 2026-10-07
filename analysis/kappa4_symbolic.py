#!/usr/bin/env python3
"""Derivacao simbolica ( sympy ) do coeficiente hiperviscoso do modo de cisalhamento no D2Q9.

Usa as recorrencias de tres momentos do Apendice A:
  - onda ao longo do eixo x:      ( m0, m1, m2 )   ->  kappa_4(0)    = a/9 [ a ( 2b - a ) - 1/4 ]
  - onda ao longo da diagonal:    ( j, sigma, r )  ->  kappa_4(pi/4) = a^2 ( b - 2a ) / 18
com a = tau - 1/2 e b = tau_nh - 1/2.

    python3 analysis/kappa4_symbolic.py
"""
import sympy as sp

a, b, q, lam = sp.symbols("a b q lambda")
tau, tnh = a + sp.Rational(1, 2), b + sp.Rational(1, 2)
C, S = sp.cos(q), sp.sin(q)


def expande(M, K2_por_q2):
    """log lambda = r2 q^2 + r4 q^4 ;  devolve ( nu, kappa_4 ) com K^2 = K2_por_q2 * q^2."""
    p = (lam * sp.eye(3) - M).det()
    r2, r4 = sp.symbols("r2 r4")
    ser = sp.series(p.subs(lam, sp.exp(r2 * q ** 2 + r4 * q ** 4)), q, 0, 5).removeO()
    s2 = sp.solve(sp.simplify(ser.coeff(q, 2)), r2)[0]
    s4 = sp.solve(sp.simplify(ser.coeff(q, 4).subs(r2, s2)), r4)[0]
    return sp.simplify(-s2 / K2_por_q2), sp.factor(s4 / K2_por_q2 ** 2)


# eixo x:  m0 = sum c_y df,  m1 = sum c_x c_y df,  m2 = sum c_x^2 c_y df
col = sp.Matrix([[1, 0, 0], [0, 1 - 1 / tau, 0], [1 / (3 * tnh), 0, 1 - 1 / tnh]])
pro = sp.Matrix([[1, -sp.I * S, C - 1], [0, C, -sp.I * S], [0, -sp.I * S, C]])
print("eixo:     nu = %s    kappa_4 = %s" % expande(pro * col, 1))

# diagonal:  j = sum ( c_y - c_x ) df,  sigma = sum ( c_x^2 - c_y^2 ) df,  r = f_6 - f_8
col = sp.Matrix([[1, 0, 0], [0, 1 - 1 / tau, 0], [1 / (6 * tnh), 0, 1 - 1 / tnh]])
pro = sp.Matrix([[C, sp.I * S, 2 - 2 * C], [sp.I * S, C, -2 * sp.I * S], [0, 0, 1]])
print("diagonal: nu = %s    kappa_4 = %s" % expande(pro * col, 2))
