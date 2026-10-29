"""多段の合体の計算（VIC、w=0.04・壁 ±2、run_w004_zw2）を見る: θ 方向の渦の数と、芯の円の平均渦度 / γ0。
w=0.08・壁 ±4 を周長だけ 2 倍にしたのと同じ相似形なので、w=0.08 の r=0.2 に当たるのは r=0.1、山を数える帯は |z|<0.3。
保存は 5 万粒子の抜き出し（等しい循環）なので、ω/γ0 は粒子の面密度 / 初期の帯の面密度。中心は r の円に入る粒子数が最大の点。

    python merge_w004.py run_w004_zw2.npz 0.04
"""
import math, sys
import numpy as np
fn, w = sys.argv[1], float(sys.argv[2])
d = np.load(fn)
ts = sorted([int(k[1:]) for k in d.files if k.startswith("t") and k[1:].isdigit()])
r = 0.2 * w / 0.08; zb = 0.6 * w / 0.08
for t in ts:
    th, z = d[f"t{t}"]
    n = th.size; rho0 = n / (2 * math.pi * 2 * w)
    # 山の数: |z|<zb の粒子の θ 分布、128 区間、9 区間の移動平均、平均の 1.3 倍を超える山
    h, _ = np.histogram(th[np.abs(z) < zb] % (2 * math.pi), bins=128, range=(0, 2 * math.pi))
    hs = np.convolve(np.r_[h[-4:], h, h[:4]], np.ones(9) / 9, mode="valid")
    pk = [i for i in range(128) if hs[i] > 1.3 * hs.mean() and hs[i] >= hs[i - 1] and hs[i] >= hs[(i + 1) % 128]]
    # 芯: 粗い格子で候補を取り、その周りで円の粒子数を最大にする点
    H, xe, ye = np.histogram2d(th % (2 * math.pi), z, bins=(256, 64), range=((0, 2 * math.pi), (-0.5, 0.5)))
    i, j = np.unravel_index(np.argmax(H), H.shape)
    c0 = (0.5 * (xe[i] + xe[i + 1]), 0.5 * (ye[j] + ye[j + 1]))
    best = (0, c0)
    for dx in np.linspace(-2 * r, 2 * r, 21):
        for dy in np.linspace(-r, r, 11):
            cx, cy = c0[0] + dx, c0[1] + dy
            dth = (th - cx + math.pi) % (2 * math.pi) - math.pi
            k = int(((dth ** 2 + (z - cy) ** 2) < r * r).sum())
            if k > best[0]:
                best = (k, (cx, cy))
    k, (cx, cy) = best
    dth = (th - cx + math.pi) % (2 * math.pi) - math.pi
    inn = lambda rr: ((dth ** 2 + (z - cy) ** 2) < rr * rr).sum() / (math.pi * rr * rr) / rho0
    print(f"t={t:4d} 山 {len(pk)} 個  芯 r={r/4:.3f}/{r/2:.3f}/{r:.3f} = {inn(r/4):.3f}/{inn(r/2):.3f}/{inn(r):.3f}  中心 θ={cx:.2f} z={cy:+.3f}")
