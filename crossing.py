"""NS のログから、芯（r=0.2）が多段 MRS(E)（eps3）を上から下へ横切る点を E で線形に内挿して出す。格子・種の比較用。

    python crossing.py ns_w008_nu3e-4.log ns_w008_nu3e-4_lores_e15.log ...
"""
import re, sys
import numpy as np

g = {"__name__": "x"}
argv, sys.argv = sys.argv, ["x"]
exec(open("fig_core_energy.py").read().split("fig, ax")[0], g)
eL = np.array([e for e, _ in g["mrsL"]]); yL = np.array([y for _, y in g["mrsL"]])
RX = r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)"
for f in argv[1:]:
    rows = [(float(m.group(1)), float(m.group(2)), float(m.group(5))) for l in open(f) for m in [re.search(RX, l)] if m and float(m.group(1)) >= 40]
    d = [(t, E, y, y - np.interp(E, eL, yL)) for t, E, y in rows if eL.min() <= E <= eL.max()]
    cr = [(a[0] + (b[0] - a[0]) * a[3] / (a[3] - b[3]), a[1] + (b[1] - a[1]) * a[3] / (a[3] - b[3]), a[2]) for a, b in zip(d, d[1:]) if a[3] > 0 >= b[3]]
    at = {t: (E, y, dd) for t, E, y, dd in d}
    s = " ".join(f"t={t:.0f} E={E:.3f} 芯={y:.3f}" for t, E, y in cr) or "なし"
    print(f"{f}: 横切り {s}  t=50 の差 {at.get(50.0, (0, 0, float('nan')))[2]:+.3f}  最後 t={rows[-1][0]:.0f} E={rows[-1][1]:.4f} 差 {d[-1][3]:+.3f}")
