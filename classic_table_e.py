"""表1（古典の帯の厚さ）の MRS を、NS のその時刻のエネルギーにそろえた値に置き換える。
MRS は E/E0=1.00・0.95・0.90・0.85 で解いた芯の円の平均（r=0.2）を E で線形に内挿する。

    python classic_table_e.py          （2値の MRS）
    python classic_table_e.py multi    （NS の初期場の渦度分布で解いた多段の MRS、mrs_levels.py）
"""
import math, re, sys
import numpy as np

R = 0.2
def circle(f, w):
    d = np.load(f); G, psi = d["G"], d["psi"]
    nth, nz1 = G.shape; ZW = 4.0
    th = np.arange(nth) * 2 * math.pi / nth; z = -ZW + np.arange(nz1) * 2 * ZW / (nz1 - 1)
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    dth = (th[:, None] - th[i] + math.pi) % (2 * math.pi) - math.pi
    return float(G[(dth**2 + (z[None, :] - z[j])**2) < R * R].mean()) * w

def ns_at(f, t):
    for line in open(f):
        m = re.match(rf"t=\s*{t}\.0 E/E0=([\d.]+) core\(r=\.05/\.1/\.2\)=[\d.]+/[\d.]+/([\d.]+)", line)
        if m:
            return float(m.group(1)), float(m.group(2))
    return None

for w in (0.107, 0.214):
    es = [1.0, 0.95, 0.90, 0.85]
    if sys.argv[1:] == ["multi"]:
        fs = [f"mrsL_w{w}_e1.0.npz"] + [f"mrsL_w{w}_e{e:.2f}.npz" for e in es[1:]]
    else:
        fs = [f"mrs_w{w}.npz"] + [f"mrs_w{w}_e{e:.2f}.npz" for e in es[1:]]
    if not all(__import__("os").path.exists(f) for f in fs):
        print(f"w={w}: MRS の解がそろっていない（{', '.join(f for f in fs if not __import__('os').path.exists(f))}）")
        continue
    cs = [circle(f, w) for f in fs]
    print(f"w={w}  MRS 芯 r=0.2: " + "  ".join(f"E={e:.2f}:{c:.3f}" for e, c in zip(es, cs)))
    for nu in ("3.3e-4", "6.7e-4", "1.3e-3"):
        row = []
        for t in (50, 100):
            r = ns_at(f"ns_w{w}_nu{nu}.log", t)
            if r is None:
                row.append("—"); continue
            E, c = r
            m = float(np.interp(E, es[::-1], cs[::-1]))
            row.append(f"NS {c:.3f} / MRS(E={E:.3f}) {m:.3f} / 差 {c - m:+.3f}")
        print(f"  ν={nu}: t=50 {row[0]}   t=100 {row[1]}")
