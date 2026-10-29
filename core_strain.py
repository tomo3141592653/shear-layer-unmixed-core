"""芯にかかる外のひずみ: 芯（ω>q·最大）を取り除いた渦度で ψ を解き、芯の重心でのひずみ e=√(((ψ_θθ−ψ_zz)/2)²+ψ_θz²) を出す。
芯の平均渦度 ω_c と縦横比 λ と並べ、e/ω_c が λ とともに増えるかを見る（Moore–Saffman 型の釣り合いの手がかり）。
芯を周期領域ごと取り除くと、芯自身の θ 方向の写し（間隔 2π の列）が作るひずみも消える。この分は写しの直接和で足す
（Guha, Rahmani & Lawrence 2013 は点渦の列で Γπ/(6λ²) と見積もった）。e は周り（ハロー）だけ、e_tot は写しを足したもの。

    python core_strain.py --file ns_w008_nu3e-4.npz
"""
import argparse, math
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--file", required=True)
ap.add_argument("--q", type=float, default=0.5)
a = ap.parse_args()
ZW = 4.0
d = np.load(a.file)
for k in sorted([k for k in d.files if k.startswith("t") and k != "t0"], key=lambda s: int(s[1:])):
    om = d[k].astype(float)[:, 1:]
    nth, M = om.shape; nz = M + 1
    dth = 2 * math.pi / nth; dz = 2 * ZW / nz
    jj = np.arange(1, M + 1)
    Sm = np.sin(math.pi * jj[:, None] * jj[None, :] / nz) * math.sqrt(2.0 / nz)
    kt = np.fft.fftfreq(nth, d=1.0 / nth); kz = math.pi * jj / (2 * ZW)
    K2 = kt[:, None] ** 2 + kz[None, :] ** 2
    m = om > a.q * om.max()
    ext = np.where(m, 0.0, om)
    psi = np.fft.ifft((np.fft.fft(ext, axis=0) @ Sm / K2) @ Sm, axis=0).real
    W = om * m
    i0 = np.argmax(om.max(1))
    th = np.arange(nth) * dth; z = -ZW + jj * dz
    dt_ = (th[:, None] - th[i0] + math.pi) % (2 * math.pi) - math.pi
    x = (W * dt_).sum() / W.sum(); y = (W * z[None, :]).sum() / W.sum()
    X = dt_ - x; Y = z[None, :] - y
    cxx = (W * X * X).sum() / W.sum(); czz = (W * Y * Y).sum() / W.sum(); cxz = (W * X * Y).sum() / W.sum()
    ev = np.linalg.eigvalsh([[cxx, cxz], [cxz, czz]]); lam = math.sqrt(ev[1] / ev[0])
    ic = int(round((th[i0] + x) % (2 * math.pi) / dth)) % nth; jc = int(np.argmin(abs(z - y)))
    # 芯の重心での 2 階微分（中心差分、θ は周期）
    p = lambda di, dj: psi[(ic + di) % nth, jc + dj]
    ptt = (p(1, 0) - 2 * p(0, 0) + p(-1, 0)) / dth ** 2
    pzz = (p(0, 1) - 2 * p(0, 0) + p(0, -1)) / dz ** 2
    ptz = (p(1, 1) - p(1, -1) - p(-1, 1) + p(-1, -1)) / (4 * dth * dz)
    e = math.hypot((ptt - pzz) / 2, ptz)
    # 芯の写し（θ 方向に 2πn、n≠0）が芯の重心に作る ψ の 2 階微分。核 −(1/2π)ln r の微分の直接和（壁の写しは e^{-8} 程度で無視）
    sel = W > 0
    Xs, Ys, Ws = np.broadcast_to(dt_ - x, W.shape)[sel], np.broadcast_to(Y, W.shape)[sel], W[sel] * dth * dz
    itt = izz = itz = 0.0
    for n in list(range(-400, 0)) + list(range(1, 401)):
        rx = -(Xs + 2 * math.pi * n); ry = -Ys; r4 = (rx * rx + ry * ry) ** 2
        itt += -(Ws * (ry * ry - rx * rx) / r4).sum() / (2 * math.pi)
        itz += (Ws * 2 * rx * ry / r4).sum() / (2 * math.pi)
    izz = -itt
    Ptt, Pzz, Ptz = ptt + itt, pzz + izz, ptz + itz
    e_img = math.hypot((itt - izz) / 2, itz); e_tot = math.hypot((Ptt - Pzz) / 2, Ptz)
    phi_tot = 0.5 * math.degrees(math.atan2(Pzz - Ptt, 2 * Ptz))
    Gc = Ws.sum()
    # 芯の範囲で重みつき平均したハローのひずみ（中心の一点との違いで、ひずみが芯の大きさの範囲で一様かを見る）
    Ftt = (np.roll(psi[:, 1:-1], -1, 0) - 2 * psi[:, 1:-1] + np.roll(psi[:, 1:-1], 1, 0)) / dth ** 2
    Fzz = (psi[:, 2:] - 2 * psi[:, 1:-1] + psi[:, :-2]) / dz ** 2
    Ftz = (np.roll(psi[:, 2:], -1, 0) - np.roll(psi[:, 2:], 1, 0) - np.roll(psi[:, :-2], -1, 0) + np.roll(psi[:, :-2], 1, 0)) / (4 * dth * dz)
    Wi = W[:, 1:-1]
    mtt, mzz, mtz = [(F * Wi).sum() / Wi.sum() for F in (Ftt, Fzz, Ftz)]
    e_avg = math.hypot((mtt - mzz) / 2, mtz)
    e_avg_tot = math.hypot((mtt + itt - mzz - izz) / 2, mtz + itz)
    # 外の流れのひずみの伸張軸の向き（θ 軸からの角度）: u_θ=ψ_z, u_z=−ψ_θ → s_θθ=ψ_θz, s_zz=−ψ_θz, s_θz=(ψ_zz−ψ_θθ)/2
    phi = 0.5 * math.degrees(math.atan2(pzz - ptt, 2 * ptz))
    ang = 0.5 * math.degrees(math.atan2(2 * cxz, cxx - czz))
    wc = float(om[m].mean())
    print(f"{k}: 芯の平均 ω_c={wc:.3f}  縦横比 λ={lam:.4f}  外のひずみ e={e:.4f}  e/ω_c={e / wc:.4f}  外の渦度(芯の位置) {-(ptt + pzz):.3f}  伸張軸 {phi:+.0f}° 芯の長軸 {ang:+.0f}°")
    print(f"     写しのひずみ {e_img:.4f}（点渦の列 Γπ/(6·4π²)={Gc * math.pi / (6 * 4 * math.pi**2):.4f}、芯の循環 {Gc:.3f}）  合計 e_tot={e_tot:.4f}  e_tot/ω_c={e_tot / wc:.4f}  伸張軸 {phi_tot:+.0f}°")
    print(f"     芯の範囲で重みつき平均: ハロー {e_avg:.4f}（中心 {e:.4f}）  写しを足して {e_avg_tot:.4f}  e/ω_c={e_avg_tot / wc:.4f}")
