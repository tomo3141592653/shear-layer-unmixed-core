"""粘性なし VIC のスナップショット（5万粒子の抜き出し）から、ψ 最大点のまわり半径 r の円の平均渦度 / γ0 を出す（ns_band.py の core() と同じ量）。
ψ 最大点は 64×65 の粗視化で解いた ψ から取り、その升目のまわりを細かく探して粒子密度の重心で詰める。

    python vic_core_circle.py
"""
import math, sys
import numpy as np
FILE = sys.argv[1] if len(sys.argv) > 1 else "run_w008.npz"   # 別の種: seeds_w008_s3.npz など（鍵は t60p0 の形）
sys.argv = [sys.argv[0]]
exec(open("ns_omegapsi.py").read().split("times =")[0])

md = np.load(FILE)
N, NS = 500_000, 50_000
w = 0.08; gam0 = 1 / w; Gv = gam0 * (2 * math.pi * 2 * w) / N
for t in (35, 48, 60, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900):
    k = f"t{t}" if f"t{t}" in md.files else f"t{t}p0"
    if k not in md.files:
        continue
    th, zz = md[k]
    G = from_md(th, zz)
    G = G / (G.sum() * dth * dz) * 4 * math.pi
    psi = np.zeros_like(G)
    psi[:, 1:-1] = np.fft.ifft((-(np.fft.fft(G[:, 1:-1], axis=0) @ S) / den) @ S, axis=0).real
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    c_th, c_z = i * dth + dth / 2, -ZW + j * dz
    for _ in range(5):          # 重心で中心を詰める（半径 0.15）
        d_th = (th - c_th + math.pi) % (2 * math.pi) - math.pi
        m = d_th**2 + (zz - c_z)**2 < 0.15**2
        c_th = c_th + d_th[m].mean(); c_z = zz[m].mean()
    d_th = (th - c_th + math.pi) % (2 * math.pi) - math.pi
    out = []
    for r in (0.05, 0.1, 0.2):
        n = np.sum(d_th**2 + (zz - c_z)**2 < r * r) * (N / NS)
        out.append(n * Gv / (math.pi * r * r) / gam0)
    print(f"t={t} core(r=.05/.1/.2)=" + "/".join(f"{x:.3f}" for x in out))
