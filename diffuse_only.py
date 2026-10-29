"""拡散だけの模型: NS の t=t0 の場を初期値に、移流を切って ∂ω/∂t = ν∇²ω だけで流し、エネルギーと芯（r=0.2 の円の平均）の軌跡を NS と比べる。
θ はフーリエ、z は壁（ω=0）で正弦展開なので厳密に解ける。保存場は物理側 256×256（z=-4..4−Δ の格子点）。

    python diffuse_only.py --nu 1e-4 --file ns_w008_nu1e-4_long.npz --t0 50
"""
import argparse, math, re
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--nu", type=float, required=True)
ap.add_argument("--file", required=True)
ap.add_argument("--log", default="")
ap.add_argument("--t0", type=int, default=50)
a = ap.parse_args()
w = 0.08; gam0 = 1 / w; ZW = 4.0
om0 = np.load(a.file)[f"t{a.t0}"].astype(float)            # (θ 256, z 256)、z_j = -4 + j·(8/256)
nth, nz = om0.shape
dth = 2 * math.pi / nth; dz = 2 * ZW / nz
# z=-4 の点（j=0）が壁。壁 ω=0 として内側 j=1..nz-1 を正弦展開（もう一方の壁は z=+4 = j=nz）
inner = om0[:, 1:]
M = nz - 1; jj = np.arange(1, M + 1)
Sm = np.sin(math.pi * jj[:, None] * jj[None, :] / nz) * math.sqrt(2.0 / nz)
kt = np.fft.fftfreq(nth, d=1.0 / nth)
kz = math.pi * jj / (2 * ZW)
K2 = kt[:, None] ** 2 + kz[None, :] ** 2
hat0 = np.fft.fft(inner, axis=0) @ Sm
th = np.arange(nth) * dth; z = -ZW + np.arange(1, nz) * dz

def state(t):
    h = hat0 * np.exp(-a.nu * K2 * (t - a.t0))
    om = np.fft.ifft(h @ Sm, axis=0).real
    psi = np.fft.ifft((h / K2) @ Sm, axis=0).real
    E = 0.5 * float((om * psi).sum()) * dth * dz
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    d = (th[:, None] - th[i] + math.pi) % (2 * math.pi) - math.pi
    core = float(om[(d**2 + (z[None, :] - z[j])**2) < 0.04].mean()) / gam0
    return E, core

E_ns = {}
if a.log:
    for l in open(a.log):
        m = re.search(r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)", l)
        if m: E_ns[int(float(m.group(1)))] = (float(m.group(2)), float(m.group(5)))
Eref, _ = state(a.t0)
e0 = E_ns.get(a.t0, (1, 0))[0]
print(f"ν={a.nu} t0={a.t0}  （拡散だけ: E/E0, 芯）  （NS: E/E0, 芯）")
for t in range(a.t0, 701, 50):
    E, c = state(t)
    ns = E_ns.get(t)
    print(f"t={t:4d}  拡散 {E/Eref*e0:.4f} {c:.3f}" + (f"   NS {ns[0]:.4f} {ns[1]:.3f}" if ns else ""))
