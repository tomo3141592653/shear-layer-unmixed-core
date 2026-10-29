"""超粘性の計算で、エンストロフィーと渦度の面積分布（Casimir の代わりに見る量）がどれだけ変わったかを、保存した場から出す。
保存した場は物理側（|z|<Z_W）を格子点1つおきに取り出したもの（平均ではない）なので、点の平均がそのまま面積平均の見積もりになる。

    python hyper_casimir.py ns_w008_hyper18_t100.npz ns_w008_hyper_g1024.npz ns_w008_hyper16.npz
"""
import sys
import numpy as np

BINS = [(-np.inf, 0.0), (0.0, 0.05), (0.05, 0.5), (0.5, 0.95), (0.95, 1.0), (1.0, np.inf)]   # ω/γ0 の区切り
for f in sys.argv[1:]:
    d = np.load(f)
    w = float(d["w"])
    keys = [k for k in ("t0", "t50", "t100") if k in d.files]
    z0 = None
    print(f)
    for k in keys:
        om = d[k].astype(np.float64) * w                    # ω/γ0（γ0=1/w）
        z = 0.5 * np.mean(om**2)
        z0 = z0 or z
        band = om > 0.05                                    # 帯の流体（ω/γ0>0.05）の面積で割った割合
        frac = [np.mean((om > a) & (om <= b)) / np.mean(band) for a, b in BINS]
        print(f"  {k:5s} Z/Z0={z / z0:.3f}  面積（帯の面積に対する割合）: " +
              " ".join(f"{'<0' if a == -np.inf else f'{a:g}'}〜{'' if b == np.inf else f'{b:g}'}:{p:.3f}" for (a, b), p in zip(BINS, frac)))
