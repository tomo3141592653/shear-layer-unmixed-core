"""円筒面上のせん断層（一様渦度の帯）を、粘性のある 2D Navier–Stokes の擬スペクトル法で解く。
vic.py（粘性なしの渦糸セル法）と同じ設定: R=1、θ∈[0,2π) 周期、z=±ZW に ψ=0 の滑る壁、帯 |z|<w に γ0=1/w、総循環 4π。

壁は鏡像で表す: z∈[-ZW, 3ZW) を周期にし、z'∈[ZW, 3ZW) には ω(z') = -ω(2ZW - z') を置く（壁 z=ZW について奇）。
奇の対称性は時間発展で保たれるので、z=±ZW で ψ=0・ω=0 の滑る壁と同じになる。
∇²ψ = -ω（vic.py と同じ符号）、dθ/dt = ∂ψ/∂z、dz/dt = -∂ψ/∂θ、∂ω/∂t + u·∇ω = ν∇²ω。

    uv run python ns_band.py --w 0.08 --nu 1e-5 --tmax 400 --out ns_w008_nu1e-5
"""
import argparse
import math
import time

import numpy as np
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--w", type=float, default=0.08)
ap.add_argument("--nu", type=float, default=1e-5)
ap.add_argument("--nu4", type=float, default=0.0)      # 超粘性 ν4 ∇⁸（粘性なしの近似。打ち切り波数で減衰率がほぼ1になる大きさ 1e-18 程度）
ap.add_argument("--nth", type=int, default=512)
ap.add_argument("--nz", type=int, default=1024)        # 鏡像を含む z の長さ 4ZW の点数
ap.add_argument("--zw", type=float, default=4.0)
ap.add_argument("--dt", type=float, default=0.005)
ap.add_argument("--tmax", type=float, default=400.0)
ap.add_argument("--noise", type=float, default=1e-3)   # 帯の中の渦度のゆらぎ（γ0 に対する割合）
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", default="ns")
ap.add_argument("--eps", type=float, default=3.0)       # 帯の縁をならす幅（格子間隔 dz の何倍か）
a = ap.parse_args()

dev = "cuda"
DT = torch.float64
torch.manual_seed(a.seed)
NTH, NZ, ZW, w = a.nth, a.nz, a.zw, a.w
LZ = 4 * ZW
th = torch.arange(NTH, dtype=DT, device=dev) * 2 * math.pi / NTH
z = -ZW + torch.arange(NZ, dtype=DT, device=dev) * LZ / NZ
dz = LZ / NZ
kt = torch.fft.fftfreq(NTH, d=1.0 / NTH).to(dev).to(DT)                 # θ の波数（整数）
kz = torch.fft.rfftfreq(NZ, d=LZ / NZ).to(dev).to(DT) * 2 * math.pi    # z の波数
KT, KZ = torch.meshgrid(kt, kz, indexing="ij")
K2 = KT**2 + KZ**2
K2i = torch.where(K2 > 0, 1.0 / K2, torch.zeros_like(K2))
mask = ((KT.abs() < NTH / 3) & (KZ.abs() < (NZ / 3) * 2 * math.pi / LZ)).to(DT)   # 2/3 則

gam0 = 1.0 / w
Z = z[None, :].expand(NTH, NZ)
eps = a.eps * dz                                           # 帯の縁を格子 a.eps 点ぶんでならす（既定 3）
band = 0.5 * (torch.tanh((w - Z.abs()) / eps) + 1)                    # 物理の帯 |z|<w
img = 0.5 * (torch.tanh((w - (Z - 2 * ZW).abs()) / eps) + 1)          # 鏡像 z=2ZW（符号反転）
om = gam0 * (band - img)
om = om + gam0 * a.noise * torch.randn(NTH, NZ, dtype=DT, device=dev) * (band - img).abs()
om = om * (4 * math.pi / float((om * (Z.abs() < ZW)).sum() * (2 * math.pi / NTH) * dz))   # 物理側の総循環を 4π に
omh = torch.fft.rfft2(om) * mask

def rhs(omh):
    psih = omh * K2i                                       # ∇²ψ = -ω → ψ̂ = ω̂/k²
    ps_t = torch.fft.irfft2(1j * KT * psih, s=(NTH, NZ))
    ps_z = torch.fft.irfft2(1j * KZ * psih, s=(NTH, NZ))
    om_t = torch.fft.irfft2(1j * KT * omh, s=(NTH, NZ))
    om_z = torch.fft.irfft2(1j * KZ * omh, s=(NTH, NZ))
    adv = ps_z * om_t - ps_t * om_z                        # u·∇ω（uθ=ψ_z, uz=-ψ_θ）
    return (-torch.fft.rfft2(adv) - (a.nu * K2 + a.nu4 * K2**4) * omh) * mask

def diag(omh):
    om = torch.fft.irfft2(omh, s=(NTH, NZ))
    psi = torch.fft.irfft2(omh * K2i, s=(NTH, NZ))
    phys = (Z.abs() < ZW)
    cell = (2 * math.pi / NTH) * dz
    E = 0.5 * float((om * psi * phys).sum()) * cell
    circ = float((om * phys).sum()) * cell
    return om, psi, E, circ

def core(om, psi, r):
    """芯の中心（物理側の ψ 最大点）のまわり半径 r の円の平均渦度 / γ0。coremeas.disk_field と同じ考え方。"""
    p = torch.where(Z.abs() < ZW, psi, torch.full_like(psi, -1e9))
    ic = np.unravel_index(int(torch.argmax(p)), p.shape)
    dtt = (th[:, None] - th[ic[0]] + math.pi) % (2 * math.pi) - math.pi
    m = dtt**2 + (Z - z[ic[1]])**2 < r * r
    return float(om[m].mean()) * w, float(th[ic[0]]), float(z[ic[1]])

_, _, E0, C0 = diag(omh)
print(f"w={w} nu={a.nu} nu4={a.nu4} grid={NTH}x{NZ} dt={a.dt} E0={E0:.6f} circ={C0:.6f}", flush=True)
snap_t = [0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900]
log, snaps = [], {}
nsteps = int(round(a.tmax / a.dt))
t0 = time.time()
for s in range(nsteps + 1):
    t = s * a.dt
    if s % int(round(10 / a.dt)) == 0:
        om, psi, E, C = diag(omh)
        c05, c1, c2 = (core(om, psi, r)[0] for r in (0.05, 0.1, 0.2))
        log.append((t, E / E0, C, c05, c1, c2, float(om.max()) * w))
        print(f"t={t:6.1f} E/E0={E/E0:.5f} core(r=.05/.1/.2)={c05:.3f}/{c1:.3f}/{c2:.3f} max={float(om.max())*w:.3f} {time.time()-t0:.0f}s", flush=True)
        if any(abs(t - x) < a.dt / 2 for x in snap_t):
            snaps[f"t{t:.0f}"] = om[:, (z.abs() < ZW)].float().cpu().numpy()[::2, ::2]
    if s == nsteps:
        break
    h = a.dt
    k1 = rhs(omh); k2 = rhs(omh + h / 2 * k1); k3 = rhs(omh + h / 2 * k2); k4 = rhs(omh + h * k3)
    omh = omh + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

np.savez_compressed(a.out + ".npz", log=np.array(log), E0=E0, nu=a.nu, w=w, **snaps)
print("done", flush=True)
