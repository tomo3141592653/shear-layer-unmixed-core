"""一定のひずみの中の孤立渦の定常の縦横比を、渦の動径の形（縁の切り立ち方）だけ変えて測る。
せん断層の芯で「測った e/ω_c が一様パッチの式の 1.1〜1.2 倍」になった原因が、芯のなめらかさと二次モーメントでの測り方にあるかを切り分ける。

  ∂ω/∂t + (u_v + u_s)·∇ω = −ν4 ∇⁴ω,  u_s = (e x, −e y)（伸張軸は x、外の渦度なし）
  u_v は周期箱（一辺 L）で解く渦自身の流れ。箱の端は吸収層で ω を消す（ひずみで伸びた糸が端で折り返さないように）

芯の量は core_strain.py と同じ: ω > q·最大 の領域で ω_c（平均）と λ（重みつき二次モーメントの固有値比）。
e は外から与えた値なので、e/ω_c と λ を一様パッチの式 e/ω = λ(λ−1)/((λ+1)(λ²+1)) と比べる。

    python strain_vortex.py --prof tanh --edge 0.005
"""
import argparse, math
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--prof", default="tanh", choices=["tanh", "gauss"])
ap.add_argument("--edge", type=float, default=0.02)     # tanh の縁の幅
ap.add_argument("--R", type=float, default=0.25)        # 芯の半径（tanh）／ガウスの幅
ap.add_argument("--w0", type=float, default=8.5)        # 芯の渦度
ap.add_argument("--e", type=float, default=0.23)        # 外のひずみ
ap.add_argument("--L", type=float, default=3.0)
ap.add_argument("--N", type=int, default=256)
ap.add_argument("--nu4", type=float, default=2e-9)
ap.add_argument("--tmax", type=float, default=60.0)
ap.add_argument("--q", type=float, default=0.5)
ap.add_argument("--ramp", type=float, default=40.0)  # ひずみを 0 から e まで上げる時間（ゆっくり上げて定常の枝に乗せる）
a = ap.parse_args()

N, L = a.N, a.L
x = (np.arange(N) - N / 2) * L / N
X, Y = np.meshgrid(x, x, indexing="ij")
k = np.fft.fftfreq(N, d=L / N) * 2 * math.pi
KX, KY = np.meshgrid(k, k, indexing="ij")
K2 = KX**2 + KY**2; K2i = np.where(K2 > 0, 1 / np.where(K2 > 0, K2, 1), 0)
deal = (np.abs(KX) < (2 / 3) * k.max()) & (np.abs(KY) < (2 / 3) * k.max())
r = np.hypot(X, Y)
if a.prof == "tanh":
    w = 0.5 * a.w0 * (1 - np.tanh((r - a.R) / a.edge))
else:
    w = a.w0 * np.exp(-(r / a.R) ** 2)
# 吸収層: 中心から 0.38L より外で ω を減衰
sponge = np.clip((np.maximum(np.abs(X), np.abs(Y)) - 0.38 * L) / (0.12 * L), 0, 1) ** 2 * 5.0
us, vs = a.e * X, -a.e * Y


def rhs(w, t):
    f = min(1.0, t / a.ramp) if a.ramp > 0 else 1.0
    wh = np.fft.fft2(w)
    ps = wh * K2i
    u = np.real(np.fft.ifft2(1j * KY * ps)); v = np.real(np.fft.ifft2(-1j * KX * ps))
    wx = np.real(np.fft.ifft2(1j * KX * wh * deal)); wy = np.real(np.fft.ifft2(1j * KY * wh * deal))
    return -(u + f * us) * wx - (v + f * vs) * wy - sponge * w


def core(w):
    m = w > a.q * w.max()
    W = w * m
    cx = (W * X).sum() / W.sum(); cy = (W * Y).sum() / W.sum()
    dx, dy = X - cx, Y - cy
    c = np.array([[(W * dx * dx).sum(), (W * dx * dy).sum()], [(W * dx * dy).sum(), (W * dy * dy).sum()]]) / W.sum()
    ev, vec = np.linalg.eigh(c)
    ang = math.degrees(math.atan2(vec[1, 1], vec[0, 1]))
    ang = (ang + 90) % 180 - 90
    return float(w[m].mean()), math.sqrt(ev[1] / ev[0]), ang, cx, cy


def patch_lam(rr):
    lo, hi = 1.0 + 1e-9, 2.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mid * (mid - 1) / ((mid + 1) * (mid * mid + 1)) < rr:
            lo = mid
        else:
            hi = mid
    return lo


dx = L / N
umax = a.w0 * a.R + a.e * L / 2
dt = 0.3 * dx / umax
kd = np.exp(-a.nu4 * K2**2 * dt)
t = 0.0; nxt = 0.0
hist = []
while t < a.tmax - 1e-9:
    k1 = rhs(w, t); k2 = rhs(w + 0.5 * dt * k1, t + 0.5 * dt); k3 = rhs(w + 0.5 * dt * k2, t + 0.5 * dt); k4 = rhs(w + dt * k3, t + dt)
    w = w + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    w = np.real(np.fft.ifft2(np.fft.fft2(w) * kd))
    t += dt
    if t >= nxt - 1e-9:
        wc, lam, ang, cx, cy = core(w)
        hist.append((t, wc, lam, ang))
        if t >= a.ramp - 2: print(f"t={t:6.2f} ω_c={wc:.3f} λ={lam:.3f} 長軸 {ang:+.0f}° e/ω_c={a.e / wc:.4f} 式の λ={patch_lam(a.e / wc):.3f} 中心 ({cx:+.3f},{cy:+.3f})", flush=True)
        nxt += 1.0
h = np.array(hist)
late = h[h[:, 0] >= a.ramp + 0.5 * (a.tmax - a.ramp)]
r_ = a.e / late[:, 1].mean()
print(f"後半平均: λ={late[:, 2].mean():.3f}±{late[:, 2].std():.3f} 長軸 {late[:, 3].mean():+.0f}° e/ω_c={r_:.4f} 式の λ={patch_lam(r_):.3f} "
      f"比（測った e/ω_c ÷ 測った λ での式の e/ω）={r_ / (late[:, 2].mean() * (late[:, 2].mean() - 1) / ((late[:, 2].mean() + 1) * (late[:, 2].mean() ** 2 + 1))):.2f}")
