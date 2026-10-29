"""NS のある時刻の場を出発点に、移流なしの Robert–Sommeria (1992) 緩和方程式で、その時刻の Casimir（渦度の段の面積）・循環・エネルギーを
保ったまま混合エントロピーを増やし、行き着いた局所最大（＝その時刻の Casimir とエネルギーでの MRS のうち、NS の状態から登れるもの）を求める。

    ∂r_k/∂t = ∇·[ D ( ∇r_k + β(t) r_k (σ_k − ω̄) ∇ψ ) ],   ω̄ = Σ σ_k r_k,  ∇²ψ = −γ0 ω̄
    β(t) = − ∫ D ∇ω̄·∇ψ / ∫ D q2 |∇ψ|²   （dE/dt = 0 にする）,  q2 = Σ r_k σ_k² − ω̄²（局所の分散）

各段の量は保存形（流束の発散）なので面積と循環は自動で保たれる。θ は周期、壁 z=±4 では流束 0（ψ=0）。格子は NS の保存場（θ 256 × z 256）の
内側の点 j=1..255。

    python relax_levels.py --snap ns_w008_nu3e-4.npz --t 300 --K 10 --steps 40000
"""
import argparse
import math

import numpy as np
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--snap", default="")
ap.add_argument("--mrsfile", default="")   # 検証用: mrs.py の解（2値）から出発し、緩和の定常点になっているかを見る
ap.add_argument("--vic", default="")    # 検証用: VIC のスナップショット（run_w008.npz）から 2 段（σ=0,1）の出発点を作る
ap.add_argument("--t", type=int, required=True)
ap.add_argument("--K", type=int, default=10)
ap.add_argument("--steps", type=int, default=40000)
ap.add_argument("--D", type=float, default=1.0)
ap.add_argument("--out", default="")
ap.add_argument("--cfl", type=float, default=0.2)
ap.add_argument("--ecorr", type=float, default=0.0)
ap.add_argument("--exactE", action="store_true")       # 離散化のレベルで dE=0 になる β を毎歩解く（このとき --ecorr は 1 歩で戻す割合）   # エネルギーのずれを β で戻す強さ（離散化でゆっくり漂う分を打ち消す）
ap.add_argument("--every", type=int, default=2000)
ap.add_argument("--sg", action="store_true")          # Scharfetter–Gummel の流束（漂流と拡散を指数のあてはめで一緒に取る。セルのペクレ数が大きくても数値拡散を足さない）
ap.add_argument("--renorm", action="store_true")      # 段0で補う代わりに、Σ_k の流束のずれを全段に風上の割合で配る
a = ap.parse_args()

dev = "cuda"; DT = torch.float64
w = 0.08; gam0 = 1 / w; ZW = 4.0
if a.mrsfile:
    Gm = np.load(a.mrsfile)["G"]                                   # (θ 256, z 257)、行 0 と 256 が壁
    om = Gm[:, 1:256]
    nth, nz = om.shape
    dth = 2 * math.pi / nth; dz = 2 * ZW / (nz + 1)
    p1 = np.clip(om / gam0, 0, 1)
    sig = np.array([0.0, 1.0]); a.K = 2
    r = np.stack([1 - p1, p1])
elif a.vic:
    th_p, zz_p = np.load(a.vic)[f"t{a.t}"]
    H, _, _ = np.histogram2d(th_p % (2 * math.pi), zz_p, bins=[256, 256], range=[[0, 2 * math.pi], [-ZW, ZW]])
    om = (H / (H.sum() * (2 * math.pi / 256) * (2 * ZW / 256)) * 4 * math.pi)[:, 1:]
    nth, nz = om.shape
    dth = 2 * math.pi / nth; dz = 2 * ZW / (nz + 1)
    p1 = np.clip(om / gam0, 0, 1)
    sig = np.array([0.0, 1.0]); a.K = 2
    r = np.stack([1 - p1, p1])                                     # 2 段の巨視状態（粒子の濃さ＝段1の割合）
else:
    om = np.load(a.snap)[f"t{a.t}"].astype(float)[:, 1:]          # 壁の行 j=0 を落とす（θ 256 × 内側 255）
    nth, nz = om.shape
    dth = 2 * math.pi / nth; dz = 2 * ZW / (nz + 1)
    s = np.clip(om / gam0, 0, None)
    edges = np.linspace(0, s.max() + 1e-12, a.K + 1)
    idx = np.clip(np.digitize(s, edges) - 1, 0, a.K - 1)
    sig = np.array([s[idx == i].mean() if np.any(idx == i) else 0.5 * (edges[i] + edges[i + 1]) for i in range(a.K)])
    r = np.stack([(idx == i).astype(float) for i in range(a.K)])    # 出発点: 各点は一つの段に属する（NS の場をそのまま段に分けたもの）


# 出発点をわずかにならす（r=0 のままだと log が発散する）。各段の領域平均と混ぜるので各段の面積・循環は変わらない
r = 0.999 * r + 0.001 * r.mean(axis=(1, 2), keepdims=True)
sig_t = torch.tensor(sig, dtype=DT, device=dev)[:, None, None]
r = torch.tensor(r, dtype=DT, device=dev)
jj = torch.arange(1, nz + 1, dtype=DT, device=dev)
Sm = torch.sin(math.pi * jj[:, None] * jj[None, :] / (nz + 1)) * math.sqrt(2.0 / (nz + 1))
kt = torch.fft.fftfreq(nth, d=1.0 / nth).to(dev).to(DT)
kz = math.pi * jj / (2 * ZW)
K2 = kt[:, None] ** 2 + kz[None, :] ** 2
Smc = Sm.to(torch.complex128)

def psi_of(wbar):
    h = torch.fft.fft(wbar * gam0, dim=0) @ Smc
    return torch.fft.ifft((h / K2) @ Smc, dim=0).real

def grad(f):
    """面の中心での勾配（θ: i+1/2、z: j+1/2）。z の壁の外は ψ=0・r は鏡で流束 0 にする側で扱う。"""
    gth = (torch.roll(f, -1, -2) - f) / dth
    gz = torch.zeros_like(f); gz[..., :-1] = (f[..., 1:] - f[..., :-1]) / dz
    return gth, gz

def avg_th(f): return 0.5 * (f + torch.roll(f, -1, -2))
def avg_z(f):
    g = torch.zeros_like(f); g[..., :-1] = 0.5 * (f[..., 1:] + f[..., :-1]); return g

def div(Fth, Fz):
    d = (Fth - torch.roll(Fth, 1, -2)) / dth
    dzf = torch.zeros_like(Fz); dzf[..., 1:] = Fz[..., 1:] - Fz[..., :-1]; dzf[..., 0] = Fz[..., 0]
    return d + dzf / dz

def energy(wbar, psi): return 0.5 * float((wbar * gam0 * psi).sum()) * dth * dz
def entropy(r): return float(-(r.clamp_min(1e-300) * r.clamp_min(1e-300).log()).sum()) * dth * dz

wbar = (sig_t * r).sum(0); psi = psi_of(wbar)
E0 = energy(wbar, psi); C0 = float(wbar.sum()) * dth * dz
dt = a.cfl * min(dth, dz) ** 2 / a.D
print(f"t={a.t} K={a.K} σ=" + " ".join(f"{x:.3f}" for x in sig) + f"  E={E0:.5f} dt={dt:.2e}", flush=True)

def core(wbar, psi, R=0.2):
    i, j = np.unravel_index(int(torch.argmax(psi)), psi.shape)
    th = torch.arange(nth, dtype=DT, device=dev) * dth; z = -ZW + torch.arange(1, nz + 1, dtype=DT, device=dev) * dz
    d = (th[:, None] - th[i] + math.pi) % (2 * math.pi) - math.pi
    return float(wbar[(d ** 2 + (z[None, :] - z[j]) ** 2) < R * R].mean())

def Bern(x):
    """B(x) = x/(e^x − 1)。x→0 で 1。"""
    small = x.abs() < 1e-6
    return torch.where(small, 1 - x / 2, x / torch.expm1(torch.where(small, torch.ones_like(x), x)))

def sg_fluxes(b, r, gp_th, gp_z, wf_th, wf_z):
    """段ごとの流束 F = −D∇r + v r を面の上で SG で取る。面 i+1/2 で F = (D/Δ)[B(−Pe) r_i − B(Pe) r_{i+1}]、Pe = vΔ/D。"""
    vth = -a.D * b * gam0 * (sig_t - wf_th[None]) * gp_th[None]
    vz = -a.D * b * gam0 * (sig_t - wf_z[None]) * gp_z[None]
    Pth = vth * dth / a.D; Pz = vz * dz / a.D
    r_r = torch.roll(r, -1, -2)
    Fth = a.D / dth * (Bern(-Pth) * r - Bern(Pth) * r_r)
    r_n = torch.zeros_like(r); r_n[..., :-1] = r[..., 1:]
    Fz = a.D / dz * (Bern(-Pz) * r - Bern(Pz) * r_n); Fz[..., -1] = 0
    wth = 0.5 * (r + r_r); wz = 0.5 * (r + r_n)                        # Σ_k の流束のずれを面の割合で配り直す
    Fth = Fth - wth * Fth.sum(0) / wth.sum(0).clamp_min(1e-300)
    Fz = Fz - wz * Fz.sum(0) / wz.sum(0).clamp_min(1e-300); Fz[..., -1] = 0
    rate = float(((Bern(-Pth) + Bern(Pth)) * a.D / dth ** 2).max() + ((Bern(-Pz) + Bern(Pz)) * a.D / dz ** 2).max())
    return Fth, Fz, rate

beta_prev = None
for step in range(a.steps + 1):
    wbar = (sig_t * r).sum(0); psi = psi_of(wbar)
    gp_th, gp_z = grad(psi)
    gw_th, gw_z = grad(wbar)
    rft, rfz = avg_th(r), avg_z(r)
    q2t = (sig_t ** 2 * rft).sum(0) - (sig_t * rft).sum(0) ** 2
    q2z = (sig_t ** 2 * rfz).sum(0) - (sig_t * rfz).sum(0) ** 2
    num = float((gw_th * gp_th).sum() + (gw_z * gp_z).sum())
    den = float((q2t * gp_th ** 2).sum() + (q2z * gp_z ** 2).sum())
    beta = -num / den / gam0 if den > 0 else 0.0
    if not a.exactE:
        beta += a.ecorr * max(abs(beta), 1.0) * (energy(wbar, psi) / E0 - 1)   # 補正は β の大きさに比例（段が多いと β が大きい）
    if step % a.every == 0:
        print(f"step={step} E/E0={energy(wbar, psi)/E0:.6f} C/C0={float(wbar.sum())*dth*dz/C0:.6f} S={entropy(r):.5f} β={beta:.3f} "
              f"芯 r=0.2 {core(wbar, psi):.3f} 最大ω/γ0 {float(wbar.max()):.3f}", flush=True)
    if a.sg:
        rft2, rfz2 = avg_th(r), avg_z(r)
        wf_th, wf_z = (sig_t * rft2).sum(0), (sig_t * rfz2).sum(0)
        dEf = lambda F: -float((psi * gam0 * (sig_t * div(F[0], F[1])).sum(0)).sum()) * dth * dz
        target = (E0 - energy(wbar, psi)) / dt * a.ecorr
        b0 = beta_prev if beta_prev is not None else beta
        b1 = b0 * 1.01 + 1e-3
        f0 = dEf(sg_fluxes(b0, r, gp_th, gp_z, wf_th, wf_z)) - target
        for _ in range(3):                                            # 割線法で離散の dE = 目標 を満たす β を解く
            f1 = dEf(sg_fluxes(b1, r, gp_th, gp_z, wf_th, wf_z)) - target
            if f1 == f0 or abs(f1) < 1e-12 * max(1.0, abs(E0)): break
            b0, b1, f0 = b1, b1 - f1 * (b1 - b0) / (f1 - f0), f1
        beta = b1; beta_prev = beta
        Fth, Fz, rate = sg_fluxes(beta, r, gp_th, gp_z, wf_th, wf_z)
        h = min(dt, a.cfl / rate)
        if step % a.every == 0:
            print(f"   SG: 解いた β={beta:.3f} 残差 dE={f1:.2e} 刻み={h:.2e}", flush=True)
        r = r - h * div(Fth, Fz)
        if step % a.every == 0:
            print(f"   最小 {float(r.min()):.2e}", flush=True)
        r = r.clamp_min(0)
        continue
    gr_th, gr_z = grad(r)
    rf_th, rf_z = avg_th(r), avg_z(r)                                 # 面の上の割合（Σ_k = 1 を保つ）
    wf_th, wf_z = (sig_t * rf_th).sum(0), (sig_t * rf_z).sum(0)       # 面の上の ω̄ を面の割合から作る → Σ_k r(σ−ω̄) = 0 が厳密
    def fluxes(b):
        vth = -a.D * b * gam0 * (sig_t - wf_th[None]) * gp_th[None]   # 段ごとの漂流の速さ（面の上）
        vz = -a.D * b * gam0 * (sig_t - wf_z[None]) * gp_z[None]
        r_up_th = torch.where(vth > 0, r, torch.roll(r, -1, -2))      # 風上の値（θ: 面 i+1/2 の左は i、右は i+1）
        r_dn = torch.zeros_like(r); r_dn[..., :-1] = r[..., 1:]
        r_up_z = torch.where(vz > 0, r, r_dn)
        Fth = vth * r_up_th; Fz = vz * r_up_z; Fz[..., -1] = 0
        if a.renorm:                                                   # Σ_k の流束のずれを風上の割合で配り直す（r=0 の段には乗らず、正の値を保ちやすい）
            Fth = Fth - r_up_th * Fth.sum(0) / r_up_th.sum(0).clamp_min(1e-300)
            Fz = Fz - r_up_z * Fz.sum(0) / r_up_z.sum(0).clamp_min(1e-300); Fz[..., -1] = 0
        fluxes.vmax = max(float(vth.abs().max()) / dth, float(vz.abs().max()) / dz)
        return Fth, Fz
    if a.exactE:
        # 流束は β の一次（風上の向きは β の符号で決まる）。離散の dE = ∫ψ γ0 dω̄ を拡散の分 A と β=±1 の漂流の分 B に分け、A + |β| B = 0 から β を解く
        Fd_th, Fd_z = -a.D * gr_th, -a.D * gr_z; Fd_z[..., -1] = 0
        dE = lambda Fth, Fz: -float((psi * gam0 * (sig_t * div(Fth, Fz)).sum(0)).sum()) * dth * dz
        A = dE(Fd_th, Fd_z) - (E0 - energy(wbar, psi)) / dt * a.ecorr   # A + βB = 目標の変化（E0 へ戻す分）を満たす β を解く
        sol = None
        for sgn in ((-1.0, 1.0) if beta < 0 else (1.0, -1.0)):
            B = dE(*fluxes(sgn))
            if B != 0 and -A / B >= 0:                                 # |β| = −A/B が正なら、その符号の風上で整合
                sol = sgn * (-A / B); break
        beta = sol if sol is not None else 0.0
        solved_beta = beta
    Fv_th, Fv_z = fluxes(beta)
    if step % a.every == 0 and a.exactE:
        print(f"   解いた β={beta:.3f} 刻み={min(dt, 0.3 / max(fluxes.vmax, 1e-30)):.2e}", flush=True)
    Fth = -a.D * gr_th + Fv_th
    Fz = -a.D * gr_z + Fv_z
    Fz[..., -1] = 0                                                 # 上の壁で流束 0（下の壁は div で 0 扱い）
    if not a.renorm:
        Fth[0] = -Fth[1:].sum(0); Fz[0] = -Fz[1:].sum(0)            # 段0は残りの補い: Σ_k の流束を厳密に 0 にして Σ_k r_k = 1 を保つ（段0が負になりうる）
    h = min(dt, 0.3 / max(fluxes.vmax, 1e-30))                     # 漂流の CFL で刻みを毎歩縮める
    r = r - h * div(Fth, Fz)
    if step % a.every == 0:
        neg = r.clamp_max(0)
        print(f"   負の量 段0 {-float(neg[0].sum())*dth*dz:.2e} 他 {-float(neg[1:].sum())*dth*dz:.2e} 最小 {float(r.min()):.2e} h={h:.2e}", flush=True)
    r = r.clamp_min(0)

if a.out:
    np.savez_compressed(a.out, r=r.cpu().numpy(), sig=sig, wbar=wbar.cpu().numpy(), psi=psi.cpu().numpy())
