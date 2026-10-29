"""多段の Miller–Robert–Sommeria 平衡。NS の初期場（帯の縁を tanh でならしたもの）が実際に持つ渦度の分布（Casimir）で解く。
2値の mrs.py は water-bag（縁が切り立った帯）の分布で解くが、NS は縁のならし eps·dz のため中間の渦度が多く（eps=3 で循環の 86%）、
2値の MRS は NS の比べる相手として正しくない（論文の査読1回目の分析(1)）。

    r_k(x) = exp(β σ_k ψ(x) + m_k) / Σ_j exp(β σ_j ψ(x) + m_j),   ω = γ0 Σ_k σ_k r_k
    ∫ r_k dA = A_k（各段の面積）、エネルギー = E（目標）

段は初期の帯の渦度 s(z) = ½(tanh((w−|z|)/ε)+1) を K 段に分けたもの（σ_k は段の中の平均、A_k は面積）。K=2・ε→0 で mrs.py の
フェルミ–ディラック形に戻る。m_k は各段の面積を一段ずつ合わせる反復（Sinkhorn 型）、β はエネルギーの二分法、ψ とは減衰つきの反復で自己無撞着に解く。

    python mrs_levels.py --eps 3 --K 12 --out mrsL_eps3
    python mrs_levels.py --eps 0.05 --K 2 --out mrsL_check   （2値の検証）
"""
import argparse
import math

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--w", type=float, default=0.08)
ap.add_argument("--nth", type=int, default=256)
ap.add_argument("--nz", type=int, default=256)
ap.add_argument("--zw", type=float, default=4.0)
ap.add_argument("--eps", type=float, default=3.0)      # NS の縁のならしの幅（NS の格子 dz=16/1024 の何倍か）
ap.add_argument("--K", type=int, default=12)
ap.add_argument("--efrac", type=float, default=1.0)
ap.add_argument("--out", default="")
a = ap.parse_args()

NTH, NZ, ZW = a.nth, a.nz, a.zw
dth = 2 * math.pi / NTH
dz = 2 * ZW / NZ
M = NZ - 1
j = np.arange(1, M + 1)
S = np.sin(math.pi * j[:, None] * j[None, :] / NZ) * math.sqrt(2.0 / NZ)
lam_z = (2 * np.cos(math.pi * j / NZ) - 2) / dz**2
kk = np.fft.fftfreq(NTH, d=1.0 / NTH)
lam_t = (2 * np.cos(2 * math.pi * kk / NTH) - 2) / dth**2
denom = lam_t[:, None] + lam_z[None, :]
TH = np.arange(NTH) * dth
Z = -ZW + np.arange(NZ + 1) * dz
dA = dth * dz
w = a.w
gam0 = 1.0 / w
TOT = 4 * math.pi


def poisson(G):
    Gh = np.fft.fft(G[:, 1:-1], axis=0) @ S
    p = np.fft.ifft((-Gh / denom) @ S, axis=0).real
    out = np.zeros((NTH, NZ + 1))
    out[:, 1:-1] = p
    return out


# 段の分布: 初期の帯の s(z) をこの 2 次元の格子に置き（エネルギーを測るのと同じ場）、内側の升目の渦度を K 段に分ける。
# 段・面積・循環・エネルギーがすべて同じ格子の上でそろう
eps = a.eps * 16.0 / 1024
zf = np.linspace(-ZW, ZW, 400001)
sf = 0.5 * (np.tanh((w - np.abs(zf)) / eps) + 1)
G0 = gam0 * np.interp(Z, zf, sf)[None, :] * np.ones((NTH, 1))
G0[:, 0] = G0[:, -1] = 0
G0 *= TOT / (G0.sum() * dA)
E0 = 0.5 * (G0 * poisson(G0)).sum() * dA * a.efrac
s_in = (G0[:, 1:-1] / gam0).ravel()
edges = np.linspace(0, s_in.max() + 1e-12, a.K + 1)
edges[1] = min(edges[1], 0.02)                                  # 段0 は背景（ごく弱い裾を含む）
idx = np.clip(np.digitize(s_in, edges) - 1, 0, a.K - 1)
sig = np.array([s_in[idx == k].mean() if np.any(idx == k) else 0.5 * (edges[k] + edges[k + 1]) for k in range(a.K)])
A = np.array([(idx == k).sum() * dA for k in range(a.K)])
keep = A > 0                                                     # 面積 0 の段（その渦度の升目が無い）は除く（log 0 を避ける）
sig, A = sig[keep], A[keep]
a.K = len(sig)
circ = gam0 * (sig * A).sum()
print(f"eps={a.eps}dz K={a.K} σ=" + " ".join(f"{x:.3f}" for x in sig) + f"  面積=" + " ".join(f"{x:.4f}" for x in A) + f"  循環/4π={circ/TOT:.5f}")
print(f"E0(grid)={E0:.5f}", flush=True)

sig_c = sig[:, None, None]
Acol = A[:, None]


def levels(psi, beta, m):
    x = beta * sig_c * psi[None] + m[:, None, None]
    x -= x.max(0, keepdims=True)
    e = np.exp(x)
    return e / e.sum(0, keepdims=True)


def fit_m(psi, beta, m, n=5):
    """各段の面積が A_k になるよう m_k を合わせる（Sinkhorn 型の交互の補正）。壁の行は数えない。"""
    for _ in range(n):
        r = levels(psi, beta, m)
        got = r[:, :, 1:-1].sum(axis=(1, 2)) * dA
        m = m + np.log(A / np.maximum(got, 1e-300))
        m -= m[0]
    return m


def entropy(r):
    rr = np.clip(r[:, :, 1:-1], 1e-300, 1)
    return float(-(rr * np.log(rr)).sum() * dA)


def relax(beta, G, m, n_it=400):
    psi = poisson(G)
    for _ in range(n_it):
        m = fit_m(psi, beta, m)
        r = levels(psi, beta, m)
        Gn = gam0 * (sig_c * r).sum(0)
        Gn[:, 0] = Gn[:, -1] = 0
        G = 0.7 * G + 0.3 * Gn
        psi = poisson(G)
    r = levels(psi, beta, m)
    return G, psi, r, 0.5 * (G * psi).sum() * dA, m


# 初期推測: mrs.py と同じ θ=π に寄せた楕円の芯
A_core = 2 * math.pi * 2 * w
ax, bz = math.sqrt(A_core / math.pi * 4), math.sqrt(A_core / math.pi / 4)
G_init = np.where(((TH[:, None] - math.pi) / ax) ** 2 + (Z[None, :] / bz) ** 2 <= 1, gam0, 0.0)
G_init *= TOT / (G_init.sum() * dA)

lo, hi = math.log(0.5), math.log(8.0)                            # これまでの解は β=2.4〜3.8。幅を狭めて 16 回（β の精度 0.02%）
m0 = np.log(A / A.sum())
for b in range(16):
    beta = math.exp(0.5 * (lo + hi))
    G, psi, r, E, m = relax(beta, G_init, m0.copy())
    got = r[:, :, 1:-1].sum(axis=(1, 2)) * dA
    print(f"beta={beta:.4f} E/E0={E/E0:.5f} 循環/4π={G.sum()*dA/TOT:.5f} 面積のずれ最大={np.max(np.abs(got/A-1)):.1e} "
          f"θ方向の振れ={np.ptp(G.sum(1))/G.sum(1).mean():.3f} 最大ω/γ0={G.max()/gam0:.3f} S={entropy(r):.5f}", flush=True)
    if E > E0:
        hi = math.log(beta)
    else:
        lo = math.log(beta)
if a.out:
    np.savez_compressed(a.out, Z=Z, G=G, psi=psi, r=r, sig=sig, A=A, beta=beta, m=m, E0=E0, S=entropy(r))
