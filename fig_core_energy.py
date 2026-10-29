"""図2: 芯の円の平均（r=0.2）をエネルギー E/E0 に対して描く。粘性ありの NS の軌跡と、同じエネルギーで解いた MRS の予言を重ねる。
比べる相手は NS の初期場の渦度分布（縁のならし eps）で解いた多段の MRS（mrs_levels.py、512 格子）。2値の MRS（256 格子）は参考に薄く描く。
両者が交わる点が、エネルギーをそろえたときの「芯が MRS を通る」点。

    python fig_core_energy.py
"""
import sys
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # 例: _e1（縁のならし eps=1 の計算）
import math, re, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = 0.2
def circle(f):
    d = np.load(f); G, psi = d["G"], d["psi"]
    nth, nz1 = G.shape; ZW = 4.0
    th = np.arange(nth) * 2 * math.pi / nth; z = -ZW + np.arange(nz1) * 2 * ZW / (nz1 - 1)
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    dth = (th[:, None] - th[i] + math.pi) % (2 * math.pi) - math.pi
    return float(G[(dth**2 + (z[None, :] - z[j])**2) < R * R].mean()) * 0.08

mrs = []
for e, f in ((1.0, "mrs_w008_256.npz"), (0.95, "mrs_w008_e0.95.npz"), (0.938, "mrs_w008_e0938.npz"), (0.90, "mrs_w008_e0.90.npz"),
             (0.85, "mrs_w008_e0.85.npz"), (0.80, "mrs_w008_e0.80.npz")):
    if os.path.exists(f):
        mrs.append((e, circle(f)))
mrs.sort()
print("2値 MRS(E):", " ".join(f"{e}:{y:.3f}" for e, y in mrs))
LE = "eps1" if SUF == "_e1" else "eps3"
mrsL = []
for e in (1.0, 0.95, 0.90, 0.85):
    f = f"mrsL_{LE}.npz" if e == 1.0 else f"mrsL_{LE}_e{e:.2f}.npz"
    if os.path.exists(f):
        mrsL.append((e, circle(f)))
mrsL.sort()
print("多段 MRS(E):", " ".join(f"{e}:{y:.3f}" for e, y in mrsL))

fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
ax.plot([e for e, _ in mrs], [y for _, y in mrs], "s:", color="#aaa", ms=3, lw=1, label="2-level MRS (water-bag)")
ax.plot([e for e, _ in mrsL], [y for _, y in mrsL], "s--", color="#333", ms=4, lw=1.2, label="multi-level MRS (smoothed band)")
for nu, f, c in (("1e-4", ("ns_w008_nu1e-4_long.log" if not SUF else f"ns_w008_nu1e-4{SUF}.log"), "#2a9d8f"), ("2e-4", f"ns_w008_nu2e-4{SUF}.log", "#8ab17d"),
                 ("3e-4", f"ns_w008_nu3e-4{SUF}.log", "#e9c46a"), ("1e-3", f"ns_w008_nu1e-3{SUF}.log", "#e76f51")):
    E, y = [], []
    for l in open(f):
        m = re.search(r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)", l)
        if m and float(m.group(1)) >= 40:
            E.append(float(m.group(2))); y.append(float(m.group(5)))
    ax.plot(E, y, "-", lw=1.5, color=c, label="NS ν=" + nu.replace("e-4", "×10⁻⁴").replace("e-3", "×10⁻³") + " (t≥40)")
    me = np.interp([e for e, _ in mrsL], E[::-1], y[::-1], left=np.nan, right=np.nan)
    cross = [(e, a, b) for (e, b), a in zip(mrsL, me) if not np.isnan(a)]
    print(nu, "NS−多段MRS at MRS energies:", " ".join(f"{e}:{a-b:+.3f}" for e, a, b in cross))
if not SUF:                                                    # ν=3e-4 を縁の物理の幅をそろえた半分・倍の格子で（2.6 節の格子の幅）
    for f, ls in (("ns_w008_nu3e-4_lores_e15.log", (0, (4, 2))), ("ns_w008_nu3e-4_hires_e6.log", (0, (1, 1.5)))):
        if not os.path.exists(f):
            continue
        E, y = [], []
        for l in open(f):
            m = re.search(r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)", l)
            if m and float(m.group(1)) >= 40:
                E.append(float(m.group(2))); y.append(float(m.group(5)))
        lab = "NS ν=3×10⁻⁴, half grid" if "lores" in f else "NS ν=3×10⁻⁴, doubled grid"
        ax.plot(E, y, color="#e9c46a", lw=1.2, ls=ls, label=lab)
ax.axhline(1.0, ls=":", color="#999", lw=1); ax.text(0.995, 1.02, "core–halo 1.00 (inviscid VIC 0.98–1.03)", fontsize=8, ha="left")
ax.set_xlabel("E / E0"); ax.set_ylabel(f"mean ω/γ0 within r={R}"); ax.set_xlim(0.8, 1.0); ax.set_ylim(0, 1.1)
ax.invert_xaxis(); ax.legend(fontsize=8, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94))
fig.tight_layout(); fig.savefig(f"fig_core_energy{SUF}.png")
