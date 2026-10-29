"""図: 粘性なしの値への連続性。芯の円の平均（r=0.2）を ν に対して、t=30・50・100 で描く。ν=0 は超粘性（同じ格子・同じ初期条件）。
小さい ν の 2 点（1e-4・2e-4）を通る直線を ν=0 まで延ばし、超粘性の値と比べる。

    python fig_nu_limit.py
"""
import re, sys
SUF = sys.argv[1] if len(sys.argv) > 1 else ""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def core(f, t):
    for line in open(f):
        m = re.match(rf"t=\s*{t}\.0 .*core\(r=\.05/\.1/\.2\)=[\d.]+/[\d.]+/([\d.]+)", line)
        if m:
            return float(m.group(1))

nus = [1e-4, 2e-4, 3e-4, 1e-3]
fig, ax = plt.subplots(figsize=(5.6, 4.0))
for t, c in ((30, "#7a9cc6"), (50, "#2f5d8a"), (100, "#14304d")):
    y0 = core(f"ns_w008_hyper16{SUF}.log", t)
    ys = [core(f"ns_w008_nu{n:.0e}.log".replace("e-04", "e-4").replace("e-03", "e-3").replace(".log", f"{SUF}.log"), t) for n in nus]
    xs = [n * 1e4 for n in nus]
    ax.plot([0] + xs, [y0] + ys, "o", color=c, ms=6, label=f"t={t}")
    s = (ys[1] - ys[0]) / (xs[1] - xs[0])
    ax.plot([0, 3.2], [ys[0] - s * xs[0], ys[0] + s * (3.2 - xs[0])], "-", color=c, lw=1, alpha=0.7)
ax.set_xlabel(r"$\nu\ (\times 10^{-4})$;  $\nu=0$: hyperviscous $\nu_4=10^{-16}$")
ax.set_ylabel(r"core mean, $r=0.2$ ($\omega/\gamma_0$)")
ax.set_xlim(-0.3, 10.5)
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"fig_nu_limit{SUF}.png", dpi=150)
for t in (30, 50, 100):
    y0 = core(f"ns_w008_hyper16{SUF}.log", t); y1 = core(f"ns_w008_nu1e-4{SUF}.log", t); y2 = core(f"ns_w008_nu2e-4{SUF}.log", t)
    print(f"t={t}: 超粘性 {y0:.3f}  1e-4・2e-4 の直線を ν=0 へ {2*y1 - y2:.3f}")
