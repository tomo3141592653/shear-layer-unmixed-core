"""論文の式(14)〜(20): コアハロー分布の自己無撞着解。vic.py と同じ格子・同じ壁で解く。

ω(ψ) = γ0            (ψ ≥ εF)
     = c (ψ − εh)     (εh < ψ < εF)
     = 0              (ψ ≤ εh)
εh = ψ(コア中心のθ, z=zmax)。総循環 4π とエネルギー E0 を保つように εF と c を更新する。

    uv run python corehalo.py --w 0.08 --zmax 2.1
"""
import argparse
import math

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--w", type=float, default=0.08)
ap.add_argument("--zmax", type=float, default=2.1)
ap.add_argument("--nth", type=int, default=256)
ap.add_argument("--nz", type=int, default=256)
ap.add_argument("--zw", type=float, default=4.0)
ap.add_argument("--iters", type=int, default=20000)
ap.add_argument("--out", default="")
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
# 初期状態のエネルギー（格子上の water-bag）。vic.py の E0 と同じ定義
G0 = np.where(np.abs(Z)[None, :] <= w, gam0, 0.0) * np.ones((NTH, 1))
# 格子が w を粗くしか分解しないので、層の総循環が 4π になるよう面積で正規化
G0 *= 4 * math.pi / (G0.sum() * dA)
E0 = 0.5 * (G0 * poisson(G0)).sum() * dA
TOT = 4 * math.pi
print(f"w={w} E0(grid)={E0:.5f}")

# 初期推測: θ=π, z=0 の楕円の芯（面積 = 層の面積の半分）＋薄い一様層
th0 = math.pi
A_core = 0.5 * 2 * math.pi * 2 * w
ax, bz = math.sqrt(A_core / math.pi * 6), math.sqrt(A_core / math.pi / 6)
G = np.where(((TH[:, None] - th0) / ax) ** 2 + (Z[None, :] / bz) ** 2 <= 1, gam0, 0.0)
G *= TOT / (G.sum() * dA)
psi = poisson(G)
iz = int(round((a.zmax + ZW) / dz))
G_init = G.copy()


def relax(fc, G, n_it=500):
    """芯に入る循環の割合 fc を固定。εF は {ψ≥εF}×γ0 の循環が fc×4π になる高さ、c は残りをハローに配る。"""
    psi = poisson(G)
    ncore = int(round(fc * TOT / (gam0 * dA)))
    for _ in range(n_it):
        ic0 = int(np.argmax(psi.max(1)))
        eh = psi[ic0, iz]
        inner = psi[:, 1:-1].ravel()
        epsF = np.partition(inner, -ncore)[-ncore] if ncore > 0 else inner.max() + 1
        core = psi >= epsF
        halo = (psi > eh) & ~core
        core[:, 0] = core[:, -1] = False; halo[:, 0] = halo[:, -1] = False
        rest = TOT - gam0 * core.sum() * dA
        hw = ((psi - eh) * halo).sum() * dA
        c = max(rest, 0.0) / hw if hw > 0 else 0.0
        Gn = np.where(core, gam0, np.where(halo, c * (psi - eh), 0.0))
        G = 0.6 * G + 0.4 * Gn
        psi = poisson(G)
    E = 0.5 * (G * psi).sum() * dA
    return G, psi, E, epsF, eh, c, rest / TOT


# fc を大きくすると芯が太り E が増える。E(fc) を二分法で E0 に合わせる。毎回、局在した初期推測から出発する
lo, hi = 0.0, 1.0
for b in range(20):
    q = 0.5 * (lo + hi)
    G, psi, E, epsF, eh, c, hfrac = relax(q, G_init)
    print(f"芯の割合={q:.5f} E/E0={E/E0:.5f} epsF={epsF:.4f} eh={eh:.4f} c={c:.4f} ハローの割合={hfrac:.3f} "
          f"θ方向の振れ={np.ptp(G.sum(1))/G.sum(1).mean():.3f}", flush=True)
    if E > E0:
        hi = q
    else:
        lo = q
Ntot = G.sum() * dA
print(f"最終 N/N0={Ntot/TOT:.6f} E/E0={E/E0:.6f}")
# 横方向の渦の割合 n(z)
nz = G.sum(0) * dth / TOT
if a.out:
    np.savez_compressed(a.out, Z=Z, nz=nz, G=G, psi=psi, epsF=epsF, eh=eh, c=c, E0=E0)
