"""円筒面上のせん断層（water-bag）の渦糸セル法シミュレーション。Pakter–Levin arXiv:2610.13615 の再現用。

単位: R=1, u0=1。θ∈[0,2π) 周期、z∈[-ZW,ZW] で ψ=0 の壁。
運動方程式（論文の式13）: dθ/dt = ∂ψ/∂z, dz/dt = -∂ψ/∂θ。∇²ψ = -Γ。

    uv run python vic.py --w 0.08 --tmax 900 --out run_w008
"""
import argparse
import math
import time

import numpy as np
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--w", type=float, default=0.08)
ap.add_argument("--np", type=int, default=500_000)
ap.add_argument("--nth", type=int, default=256)
ap.add_argument("--nz", type=int, default=256)
ap.add_argument("--zw", type=float, default=4.0)
ap.add_argument("--dt", type=float, default=0.01)
ap.add_argument("--tmax", type=float, default=900.0)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", default="run")
a = ap.parse_args()

dev = "cuda"
torch.manual_seed(a.seed)
DT = torch.float64
NTH, NZ, ZW = a.nth, a.nz, a.zw
dth = 2 * math.pi / NTH
dz = 2 * ZW / NZ                     # 格子点 z_j = -ZW + j dz, j=0..NZ（両端が壁）
M = NZ - 1                           # 内点の数

# 内点の z 方向2階差分（ディリクレ）の固有分解: 正弦変換
j = torch.arange(1, M + 1, dtype=DT, device=dev)
S = torch.sin(math.pi * j[:, None] * j[None, :] / NZ) * math.sqrt(2.0 / NZ)   # 直交行列
lam_z = (2 * torch.cos(math.pi * j / NZ) - 2) / dz**2
k = torch.fft.fftfreq(NTH, d=1.0 / NTH).to(dev).to(DT)
lam_t = (2 * torch.cos(2 * math.pi * k / NTH) - 2) / dth**2
denom = lam_t[:, None] + lam_z[None, :]           # 負

def poisson(G):
    """G: (NTH, M) 内点の渦度 → ψ (NTH, NZ+1)（壁を含む）。"""
    Gh = torch.fft.fft(G.to(torch.complex128), dim=0)
    Gh = Gh @ S.to(torch.complex128)
    ph = -Gh / denom
    ph = ph @ S.to(torch.complex128)
    p = torch.fft.ifft(ph, dim=0).real
    out = torch.zeros(NTH, NZ + 1, dtype=DT, device=dev)
    out[:, 1:-1] = p
    return out

w = a.w
gam0 = 1.0 / w
N = a.np
Gv = gam0 * (2 * math.pi * 2 * w) / N          # 1粒子の循環。総和 4π
th = torch.rand(N, dtype=DT, device=dev) * 2 * math.pi
zz = (torch.rand(N, dtype=DT, device=dev) * 2 - 1) * w

def deposit(th, zz):
    x = th / dth
    y = (zz + ZW) / dz
    i0 = torch.floor(x).long(); fx = x - i0
    j0 = torch.floor(y).long(); fy = y - j0
    i0 %= NTH; i1 = (i0 + 1) % NTH
    j0 = j0.clamp(0, NZ - 1); j1 = j0 + 1
    g = torch.zeros(NTH * (NZ + 1), dtype=DT, device=dev)
    for ii, wx in ((i0, 1 - fx), (i1, fx)):
        for jj, wy in ((j0, 1 - fy), (j1, fy)):
            g.index_add_(0, ii * (NZ + 1) + jj, wx * wy)
    return g.view(NTH, NZ + 1) * (Gv / (dth * dz))

def interp(F, th, zz):
    x = th / dth
    y = (zz + ZW) / dz
    i0 = torch.floor(x).long(); fx = x - i0
    j0 = torch.floor(y).long(); fy = y - j0
    i0 %= NTH; i1 = (i0 + 1) % NTH
    j0 = j0.clamp(0, NZ - 1); j1 = j0 + 1
    return (F[i0, j0] * (1 - fx) * (1 - fy) + F[i1, j0] * fx * (1 - fy)
            + F[i0, j1] * (1 - fx) * fy + F[i1, j1] * fx * fy)

def field(th, zz):
    G = deposit(th, zz)
    psi = poisson(G[:, 1:-1])
    psi_z = torch.zeros_like(psi)
    psi_z[:, 1:-1] = (psi[:, 2:] - psi[:, :-2]) / (2 * dz)
    psi_t = (torch.roll(psi, -1, 0) - torch.roll(psi, 1, 0)) / (2 * dth)
    return G, psi, psi_z, psi_t

def vel(th, zz):
    _, psi, pz, pt = field(th, zz)
    return interp(pz, th, zz), -interp(pt, th, zz)

def energy(th, zz):
    G, psi, _, _ = field(th, zz)
    return 0.5 * float((G * psi).sum()) * dth * dz

E0 = energy(th, zz)
print(f"w={w} N={N} E0={E0:.6f}", flush=True)
snaps_t = [0, 6, 18, 35, 48, 100, 200, 300, 400, 500, 600, 700, 800, 900]
late = [a.tmax - 10 * i for i in range(10)]     # 最後の10時刻で n(z) を平均
hist_edges = np.linspace(-3, 3, 241)
hists = []
snaps = {}
nsteps = int(round(a.tmax / a.dt))
t0 = time.time()
for s in range(nsteps + 1):
    t = s * a.dt
    if any(abs(t - x) < a.dt / 2 for x in snaps_t):
        idx = torch.randperm(N, device=dev)[:50000]
        snaps[f"t{t:.0f}"] = torch.stack([th[idx], zz[idx]]).cpu().numpy()
        print(f"t={t:.1f} E/E0-1={energy(th, zz)/E0-1:+.2e} max|z|={float(zz.abs().max()):.3f} "
              f"{time.time()-t0:.0f}s", flush=True)
    if any(abs(t - x) < a.dt / 2 for x in late):
        h, _ = np.histogram(zz.cpu().numpy(), bins=hist_edges)
        hists.append(h / N / (hist_edges[1] - hist_edges[0]))
    if s == nsteps:
        break
    h = a.dt
    k1 = vel(th, zz)
    k2 = vel(th + h / 2 * k1[0], zz + h / 2 * k1[1])
    k3 = vel(th + h / 2 * k2[0], zz + h / 2 * k2[1])
    k4 = vel(th + h * k3[0], zz + h * k3[1])
    th = (th + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])) % (2 * math.pi)
    zz = zz + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])

np.savez_compressed(a.out + ".npz", E0=E0, edges=hist_edges, nz=np.mean(hists, 0),
                    zmax_abs=float(zz.abs().max()),
                    z_q999=float(torch.quantile(zz.abs()[:1_000_000].float(), 0.999)),
                    **snaps)
print("done", flush=True)
