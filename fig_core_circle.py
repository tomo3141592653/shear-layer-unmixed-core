"""図1（改）: 芯の円の平均（ψ 最大点のまわり半径 r=0.2、γ0 単位）の時間変化。NS はログ（10 刻み）、VIC は vic_core_circle.py の値、
予言は理論の格子で同じ円の平均をとる。

    python fig_core_circle.py
"""
import sys
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # 例: _e1（縁のならし eps=1 の計算）
import math, re, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = 0.2
def log_series(f, col=4):
    t, y = [], []
    for l in open(f):
        m = re.search(r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)", l)
        if m:
            g = list(map(float, m.groups())); t.append(g[0]); y.append(g[col])
    return t, y

def theory_circle(f):
    d = np.load(f); G = d["G"]; psi = d["psi"] if "psi" in d.files else None
    nth, nz1 = G.shape; ZW = 4.0
    th = np.arange(nth) * 2 * math.pi / nth; z = -ZW + np.arange(nz1) * 2 * ZW / (nz1 - 1)
    gam0 = G.max() if psi is None else None
    i, j = np.unravel_index(np.argmax(psi if psi is not None else G), G.shape)
    dth = (th[:, None] - th[i] + math.pi) % (2 * math.pi) - math.pi
    m = dth**2 + (z[None, :] - z[j])**2 < R * R
    return float(G[m].mean()) * 0.08

vic = subprocess.run([sys.executable, "vic_core_circle.py"], capture_output=True, text=True).stdout
vt, vy = [], []
for l in vic.splitlines():
    m = re.search(r"t=(\d+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)", l)
    if m: vt.append(int(m.group(1))); vy.append(float(m.group(4)))

fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
ax.plot(vt, vy, "o-", ms=4, lw=1.5, color="#1f3b73", label="inviscid (VIC)")
for nu, f, c in (("1e-4", ("ns_w008_nu1e-4_long.log" if not SUF else f"ns_w008_nu1e-4{SUF}.log"), "#2a9d8f"), ("2e-4", f"ns_w008_nu2e-4{SUF}.log", "#8ab17d"),
                 ("3e-4", f"ns_w008_nu3e-4{SUF}.log", "#e9c46a"), ("1e-3", f"ns_w008_nu1e-3{SUF}.log", "#e76f51")):
    t, y = log_series(f)
    ax.plot(t, y, "-", lw=1.5, color=c, label="NS ν=" + nu.replace("e-4", "×10⁻⁴").replace("e-3", "×10⁻³"))
import os
hf = f"ns_w008_hyper16{SUF}.log"
if os.path.exists(hf):                                          # 同じ初期場の粘性なし（超粘性 ν4=1e-16）
    t, y = log_series(hf)
    ax.plot(t, y, "-", lw=1.5, color="#6a4c93", label="inviscid (hyperviscous, same initial field as NS)")
mlf = "mrsL_eps1.npz" if SUF == "_e1" else "mrsL_eps3.npz"       # NS の実際の初期場の渦度分布で解いた多段 MRS
for name, f, c in (("core–halo", "ch_w008.npz", "#999"), ("MRS 2-level", "mrs_w008_256.npz", "#555"), ("multi-level MRS (smoothed band)", mlf, "#b5651d")):
    y = theory_circle(f)
    ax.axhline(y, ls="--", lw=1, color=c)
    below = name in ("core–halo", "MRS 2-level")   # core–halo は上に VIC の点、2値は上に多段の線があるので下に置く
    ax.text(995, y - 0.015 if below else y + 0.012, f"{name} {y:.2f}", va="top" if below else "bottom", ha="right", fontsize=8, color=c)
    print(name, round(y, 3))
ax.set_xlabel("t"); ax.set_ylabel(f"mean ω/γ0 within r={R} of the core centre")
ax.set_xlim(0, 1000); ax.set_ylim(0, 1.1); ax.legend(fontsize=8, frameon=False, loc="lower left")
fig.tight_layout(); fig.savefig(f"fig_core_circle{SUF}.png")
