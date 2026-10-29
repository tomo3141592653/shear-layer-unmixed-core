"""図: 芯の縦横比 λ と、外のひずみと芯の渦度の比 e/ω_c。一様な楕円パッチの定常の釣り合い e/ω = λ(λ−1)/((λ+1)(λ²+1)) を重ねる。
点は core_strain.py の出力（t≥100、巻き上がり直後の揺れを除く）。e は芯自身の θ 方向の写しのひずみを含めた e_tot。

    python fig_core_strain.py
"""
import re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def pts(f):
    out = subprocess.run([sys.executable, "core_strain.py", "--file", f], capture_output=True, text=True).stdout
    r, cur = [], None                                     # 1 時刻 3 行（芯・写し・芯の範囲の平均）。行の数に頼らず t と λ を覚えて e_tot の行で拾う
    for l in out.splitlines():
        m = re.match(r"t(\d+): .*λ=([\d.]+)", l)
        if m:
            cur = (int(m.group(1)), float(m.group(2)))
        m = re.search(r"e_tot/ω_c=([\d.]+)", l)
        if m and cur and cur[0] >= 100:
            r.append((cur[1], float(m.group(1))))
    return np.array(r)

fig, ax = plt.subplots(figsize=(5.2, 4.0))
lam = np.linspace(1.0, 1.5, 200)
ax.plot(lam, lam * (lam - 1) / ((lam + 1) * (lam**2 + 1)), "-", color="#555", lw=1.2, label="uniform elliptical patch (steady)")
for f, c, lab in (("ns_w008_nu1e-4_long.npz", "#2a9d8f", r"$\nu=10^{-4}$, $t=100$–$700$"),
                  ("ns_w008_nu3e-4.npz", "#e9c46a", r"$\nu=3\times10^{-4}$, $t=100$–$400$")):
    p = pts(f)
    ax.plot(p[:, 0], p[:, 1], "o-", color=c, ms=5, lw=1, label=lab)
ax.set_xlabel(r"aspect ratio $\lambda$ of the core ($\omega>\omega_{\max}/2$)")
ax.set_ylabel(r"external strain / core vorticity  $e/\omega_c$")
ax.set_xlim(1.0, 1.5); ax.set_ylim(0, 0.1)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("fig_core_strain.png", dpi=150)
