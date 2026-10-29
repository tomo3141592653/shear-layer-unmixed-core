"""vic.py と同じ計算で、スナップショットの粒子を毎回同じ番号（最初の5万個）に固定したもの。粒子の行き先を追うため（芯の外の堀の確かめ）。

    uv run python vic_track.py --w 0.08 --tmax 300 --out track_w008
    uv run python vic_track.py --w 0.08 --tmax 300 --fine 2 --out trackf_w008   （t=100〜300 を2刻みで。鍵は t100p0 の形）
    uv run python vic_track.py --w 0.08 --tmax 210 --fine 0.1 --fine-from 200 --out trackr_w008   （芯の回転を捉える0.1刻み）

元の説明: 円筒面上のせん断層（water-bag）の渦糸セル法シミュレーション。Pakter–Levin arXiv:2610.13615 の再現用。

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
ap.add_argument("--fine", type=float, default=0, help="0 以外なら t=fine_from〜tmax をこの刻みでも記録する")
ap.add_argument("--fine-from", type=float, default=100.0)
ap.add_argument("--init-core", default="", help="a,b を与えると、層の代わりに半径 a（θ方向）・b（z方向）の一様な楕円の芯（渦度 γ0）を θ=π, z=0 に置き、残りの粒子を |z|<halo-z の帯に一様に置く（芯の外だけ）。空きの無い状態から始める実験用")
ap.add_argument("--halo-z", type=float, default=0.75)
ap.add_argument("--halo-from", default="", help="--init-core と一緒に使う。'ファイル名:枚の鍵'（例 strobe_w020:t24p0）を与えると、ハローを一様な帯でなく、その枚に記録された粒子のうち芯（ω/γ0>0.5 の升目）の外にいるものの配置から重複を許して引き、少しずらして置く。芯の中心が θ=π, z=0 に来るように平行移動する。ハローのむらを本番と同じにする実験用")
ap.add_argument("--gap-keep", type=float, default=0.0, help="--init-gap と一緒に使う。空きの範囲に落ちたハローの粒子をこの割合だけ残す（0 で空き、0.5 で半分埋まった空きから始める）")
ap.add_argument("--init-gap", action="store_true", help="--init-core と一緒に使う。芯の短軸（z 方向）の両脇、芯の大きさの 1.05〜1.6 倍・短軸から ±35° にハローを置かない（空きを入れて始める）")
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
if a.init_core:
    ca, cb = (float(v) for v in a.init_core.split(","))
    ncore = int(round(gam0 * math.pi * ca * cb / Gv))          # 芯の渦度を γ0 にする粒子数
    if ncore >= N:
        raise SystemExit("芯が大きすぎる")
    r = torch.sqrt(torch.rand(ncore, dtype=DT, device=dev)); ph = torch.rand(ncore, dtype=DT, device=dev) * 2 * math.pi
    th_c = math.pi + ca * r * torch.cos(ph); zz_c = cb * r * torch.sin(ph)
    nh = N - ncore; got = []; tot = 0
    if a.halo_from:
        import numpy as _np
        _fn, _key = a.halo_from.split(":")
        _th, _zz = _np.load(_fn + ".npz")[_key]; _th = _th % (2 * math.pi)
        _n, _dt, _dz = 64, 2 * math.pi / 64, 2 * ZW / 64
        _H, _, _ = _np.histogram2d(_th, _zz, bins=[_n, _n + 1], range=[[0, 2 * math.pi], [-ZW - _dz / 2, ZW + _dz / 2]])
        _G = _H / (_H.sum() * _dt * _dz) * 4 * math.pi * w                        # ω/γ0
        _ii = _np.clip((_th / _dt).astype(int), 0, _n - 1); _jj = _np.clip(((_zz + ZW) / _dz + 0.5).astype(int), 0, _n)
        _core = _G[_ii, _jj] > 0.5
        _tc = _np.angle(_np.exp(1j * _th[_core]).mean()) % (2 * math.pi); _zc = _zz[_core].mean()
        pool_t = torch.tensor((_th[~_core] - _tc + math.pi) % (2 * math.pi), dtype=DT, device=dev)
        pool_z = torch.tensor(_zz[~_core] - _zc, dtype=DT, device=dev)
        print(f"halo-from {a.halo_from}: 記録 {len(_th)} 個のうち芯の外 {int((~_core).sum())} 個を元にする（芯の中心 θ={_tc:.3f} z={_zc:.3f}）")
    while tot < nh:                                              # 帯に一様に置き、芯の内側に落ちたものは捨てる
        if a.halo_from:
            idx = torch.randint(0, len(pool_t), (2 * nh,), device=dev)
            t_ = (pool_t[idx] + torch.randn(2 * nh, dtype=DT, device=dev) * 0.02) % (2 * math.pi)
            z_ = pool_z[idx] + torch.randn(2 * nh, dtype=DT, device=dev) * 0.02
        else:
            t_ = torch.rand(2 * nh, dtype=DT, device=dev) * 2 * math.pi
            z_ = (torch.rand(2 * nh, dtype=DT, device=dev) * 2 - 1) * a.halo_z
        keep = ((t_ - math.pi) / ca) ** 2 + (z_ / cb) ** 2 > 1
        if a.init_gap:
            rr_ = torch.sqrt(((t_ - math.pi) / ca) ** 2 + (z_ / cb) ** 2)        # 楕円で測った半径（芯の縁が 1）
            ang_ = torch.atan2(z_ / cb, (t_ - math.pi) / ca).abs()
            gap = (rr_ > 1.05) & (rr_ < 1.6) & ((ang_ - math.pi / 2).abs() < math.radians(35))
            if a.gap_keep > 0:
                gap = gap & (torch.rand(len(gap), dtype=DT, device=dev) >= a.gap_keep)
            keep = keep & ~gap
        got.append(torch.stack([t_[keep], z_[keep]])); tot += int(keep.sum())
    h = torch.cat(got, 1)[:, :nh]
    th = torch.cat([th_c, h[0]]); zz = torch.cat([zz_c, h[1]])
    perm = torch.randperm(N, device=dev); th = th[perm]; zz = zz[perm]   # 記録する先頭5万個に芯とハローが混ざるように
    print(f"init: 芯 a={ca} b={cb} 粒子 {ncore}（循環 {ncore*Gv:.2f}/4π）、ハロー |z|<{a.halo_z} の渦度 {nh*Gv/(2*math.pi*2*a.halo_z - math.pi*ca*cb):.3f}（γ0 の {nh*Gv/(2*math.pi*2*a.halo_z - math.pi*ca*cb)/gam0:.3f} 倍）")

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
snaps_t = [0, 6, 12, 18, 24, 35, 48, 75, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900]
if a.fine:
    snaps_t = sorted(set(snaps_t) | {round(a.fine_from + a.fine * i, 6) for i in range(int(round((a.tmax - a.fine_from) / a.fine)) + 1)})
late = [a.tmax - 10 * i for i in range(10)]     # 最後の10時刻で n(z) を平均
hist_edges = np.linspace(-3, 3, 241)
hists = []
snaps = {}
nsteps = int(round(a.tmax / a.dt))
t0 = time.time()
for s in range(nsteps + 1):
    t = s * a.dt
    if any(abs(t - x) < a.dt / 2 for x in snaps_t):
        idx = torch.arange(50000, device=dev)
        snaps[f"t{t:.1f}".replace(".", "p")] = torch.stack([th[idx], zz[idx]]).cpu().numpy()
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
