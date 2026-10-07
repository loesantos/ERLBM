#!/usr/bin/env python3
"""Analise linear do D2Q9 ( BGK, RLBM e E-RLBM ) em torno do repouso.

Monta a matriz linearizada de colisao + propagacao M(K) = S(K) C para um vetor de onda K e
  1) ajusta o autovalor do modo de cisalhamento,  -log lambda = nu K^2 - kappa_4 K^4 + O(K^6),
     e compara kappa_4 com as Eqs. (15)-(17) do artigo para varias direcoes;
  2) evolui linearmente o modo de Taylor-Green ( vetores de onda ( +-k, +-k ) ) a partir de f_eq ou de
     f_eq + f^(1) e devolve ( nu_ef - nu ) / nu em 2 nu k^2 t = 1, o mesmo que bin/taylor_green --tau mede.

    python3 analysis/linear_d2q9.py
"""
import numpy as np

c = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]], float)
w = np.array([4 / 9] + [1 / 9] * 4 + [1 / 36] * 4)
cs2 = 1 / 3


def colisao(tau, tnh):
    """Matriz de colisao linearizada:  f* = E f + (1 - 1/tau) P2 (I - E) f + (1 - 1/tnh) (I - P2) (I - E) f ."""
    E = np.array([[w[i] * (1 + 3 * c[i] @ c[j]) for j in range(9)] for i in range(9)])
    H = np.einsum("ia,ib->iab", c, c) - cs2 * np.eye(2)
    P2 = np.array([[w[i] / (2 * cs2 ** 2) * np.sum(H[i] * np.outer(c[j], c[j])) for j in range(9)] for i in range(9)])
    I = np.eye(9)
    N = I - E
    return E + (1 - 1 / tau) * P2 @ N + (1 - 1 / tnh) * (I - P2) @ N


def taxa_cisalhamento(tau, tnh, K):
    C = colisao(tau, tnh)
    lam, V = np.linalg.eig(np.diag(np.exp(-1j * c @ K)) @ C)
    t = np.array([-K[1], K[0]]) / np.linalg.norm(K)
    ov = np.array([abs(t @ (c.T @ V[:, m])) / np.linalg.norm(V[:, m]) for m in range(9)])
    m = np.argmax(ov * (abs(lam) > 0.5) * abs(lam) ** 50)
    return -np.log(lam[m]).real


def ajusta_kappa4(tau, tnh, theta):
    d = np.array([np.cos(theta), np.sin(theta)])
    ks = np.linspace(0.02, 0.12, 11)
    r = np.array([taxa_cisalhamento(tau, tnh, k * d) for k in ks])
    co = np.linalg.lstsq(np.vstack([ks ** 2, ks ** 4, ks ** 6]).T, r, rcond=None)[0]
    return co[0], -co[1]


def kappa4_teorico(tau, tnh, theta):
    a, b = tau - 0.5, tnh - 0.5
    k_bar = a / 36 * (5 * a * b - 4 * a ** 2 - 0.5)
    dk = a / 72 * (6 * a * b - 1)
    return k_bar + dk * np.cos(4 * theta)


def taylor_green_linear(L, tau, tnh, init="neq"):
    k = 2 * np.pi / L
    K = np.array([k, k])
    nu = (tau - 0.5) / 3
    n = int(round(1 / (2 * nu * k * k)))
    M = np.diag(np.exp(-1j * c @ K)) @ colisao(tau, tnh)
    t = np.array([-1, 1]) / np.sqrt(2)
    f0 = w * 3 * (c @ t)
    if init == "neq":
        H = np.einsum("ia,ib->iab", c, c) - cs2 * np.eye(2)
        G = 1j * np.outer(K, t)
        f0 = f0 - tau * w / (2 * cs2 ** 2) * np.einsum("iab,ab->i", H, cs2 * (G + G.T))
    f = np.linalg.matrix_power(M, n) @ f0
    A = (t @ (c.T @ f)).real
    return -np.log(A) / (2 * k * k * n) / nu - 1


if __name__ == "__main__":
    print("kappa_4( theta ): ajuste dos autovalores x Eq. (17)")
    for tau, tnh in [(0.62, 0.62), (0.62, 1.0), (0.62, 0.65), (0.8, 0.65), (0.55, 0.65)]:
        for th in [0, np.pi / 8, np.pi / 6, np.pi / 4]:
            nu, k4 = ajusta_kappa4(tau, tnh, th)
            print("  tau=%.3f tau_nh=%.3f theta=%.3f   nu=%.6f (a/3=%.6f)   kappa4=% .6e   teoria=% .6e"
                  % (tau, tnh, th, nu, (tau - 0.5) / 3, k4, kappa4_teorico(tau, tnh, th)))
    print("\nTaylor-Green linear, L = 64: ( nu_ef - nu ) / nu   [ init f_eq | f_eq + f^(1) | a(2a-b)K^2/6 ]")
    K2 = 8 * np.pi ** 2 / 64 ** 2
    for tau in [0.55, 0.62, 0.8, 1.0]:
        for nome, tnh in [("BGK", tau), ("RLBM", 1.0), ("E-RLBM", 0.65)]:
            a, b = tau - 0.5, tnh - 0.5
            print("  tau=%.3f %-7s % .5e  % .5e  % .5e" % (tau, nome, taylor_green_linear(64, tau, tnh, "eq"),
                  taylor_green_linear(64, tau, tnh, "neq"), a * (2 * a - b) * K2 / 6))
