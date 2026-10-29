"""同じ初期条件（2値 water-bag）に対する Miller–Robert–Sommeria 平衡。
ω = γ0 / (1 + exp(-β(ψ - μ)))。総循環 4π とエネルギー E0 を保つように β, μ を決める。
格子・壁は vic.py / corehalo.py と同じ。

    uv run python mrs.py --w 0.08 --out mrs_w008
"""
import argparse
import math

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--w", type=float, default=0.08)
ap.add_argument("--nth", type=int, default=256)
ap.add_argument("--nz", type=int, default=256)
ap.add_argument("--zw", type=float, default=4.0)
ap.add_argument("--iters", type=int, default=6000)
ap.add_argument("--out", default="")
ap.add_argument("--efrac", type=float, default=1.0)   # 目標のエネルギーを初期の何倍にするか（粘性で減った場合の比較用）
ap.add_argument("--init", default="core", choices=["core", "zonal", "two"])   # 解く枝: 芯1つ／θに一様／芯2つ
a = ap.parse_args()

NTH, NZ, ZW = a.nth, a.nz, a.zw
dth = 2 * math.pi / NTH
dz = 2 * ZW / NZ
M = NZ - 1
j = np.arange(1, M + 1)
S = np.sin(math.pi * j[:, None] * j[None, :] / NZ) * math.sqrt(2.0 / NZ)
lam_z = (2 * np.cos(math.pi * j / NZ) - 2) / dz**2
k = np.fft.fftfreq(NTH, d=1.0 / NTH)
lam_t = (2 * np.cos(2 * math.pi * k / NTH) - 2) / dth**2
denom = lam_t[:, None] + lam_z[None, :]
TH = np.arange(NTH) * dth
Z = -ZW + np.arange(NZ + 1) * dz


def poisson(G):
    Gh = np.fft.fft(G[:, 1:-1], axis=0) @ S
    p = np.fft.ifft((-Gh / denom) @ S, axis=0).real
    out = np.zeros((NTH, NZ + 1))
    out[:, 1:-1] = p
    return out


w = a.w
gam0 = 1.0 / w
dA = dth * dz
TOT = 4 * math.pi
G0 = np.where(np.abs(Z)[None, :] <= w, gam0, 0.0) * np.ones((NTH, 1))
G0 *= TOT / (G0.sum() * dA)
E0 = 0.5 * (G0 * poisson(G0)).sum() * dA * a.efrac
print(f"w={w} E0(grid)={E0:.5f}")


def fd(psi, beta, mu):
    x = np.clip(-beta * (psi - mu), -60, 60)
    return gam0 / (1 + np.exp(x))


def solve_mu(psi, beta):
    """総循環が 4π になる μ を二分法で。"""
    lo, hi = psi.min() - 50 / beta, psi.max() + 50 / beta
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        n = fd(psi, beta, mid)[:, 1:-1].sum() * dA
        lo, hi = (mid, hi) if n > TOT else (lo, mid)
    return 0.5 * (lo + hi)


# 初期推測: θ=π に寄せた芯（--init two なら θ=π/2, 3π/2 に半分ずつ、zonal なら帯そのもの。対称性は緩和で保たれる）
A_core = 2 * math.pi * 2 * w
if a.init == "zonal":
    G = G0.copy()
else:
    nc = 2 if a.init == "two" else 1
    ax, bz = math.sqrt(A_core / nc / math.pi * 4), math.sqrt(A_core / nc / math.pi / 4)
    G = np.zeros((NTH, NZ + 1))
    for c in ([math.pi / 2, 3 * math.pi / 2] if nc == 2 else [math.pi]):
        G += np.where(((TH[:, None] - c) / ax) ** 2 + (Z[None, :] / bz) ** 2 <= 1, gam0, 0.0)
G *= TOT / (G.sum() * dA)
G_init = G.copy()


def entropy(G):
    """混合エントロピー −∫[p log p + (1−p) log(1−p)] dA、p = G/γ0（2値 water-bag）。"""
    p = np.clip(G / gam0, 1e-300, 1 - 1e-16)
    return float(-(p * np.log(p) + (1 - p) * np.log1p(-p)).sum() * dA)


def relax(beta, G, n_it=400):
    psi = poisson(G)
    for _ in range(n_it):
        mu = solve_mu(psi, beta)
        Gn = fd(psi, beta, mu)
        Gn[:, 0] = Gn[:, -1] = 0
        G = 0.7 * G + 0.3 * Gn
        psi = poisson(G)
    return G, psi, 0.5 * (G * psi).sum() * dA, mu


# β を大きくすると芯に集まり E が増える。log β で二分法
lo, hi = math.log(0.05), math.log(500.0)
G = G_init
for b in range(24):
    beta = math.exp(0.5 * (lo + hi))
    G, psi, E, mu = relax(beta, G_init)
    print(f"beta={beta:.4f} mu={mu:.4f} E/E0={E/E0:.5f} N/N0={G.sum()*dA/TOT:.5f} "
          f"θ方向の振れ={np.ptp(G.sum(1))/G.sum(1).mean():.3f} 最大ω/γ0={G.max()/gam0:.3f} S={entropy(G):.5f}", flush=True)
    if E > E0:
        hi = math.log(beta)
    else:
        lo = math.log(beta)
nz = G.sum(0) * dth / TOT
if a.out:
    np.savez_compressed(a.out, Z=Z, nz=nz, G=G, psi=psi, beta=beta, mu=mu, E0=E0, S=entropy(G))
