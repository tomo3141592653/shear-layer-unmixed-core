"""ε を物理的に固定した収束試験の比較: 元の格子（512×1024、ε=3dz）と倍の格子（1024×2048、ε=6dz）の芯と E/E0 を同じ時刻で並べる。
参考に、ε も半分になった倍格子（ε=3dz、ns_w008_nu3e-4_hires）と、種を変えた元の格子（ns_w008_nu3e-4_s2）の差も並べる。

    python compare_hires.py
"""
import re

RX = r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)"


def load(f):
    d = {}
    for l in open(f):
        m = re.search(RX, l)
        if m:
            d[float(m.group(1))] = (float(m.group(2)), float(m.group(5)))
    return d


base = load("ns_w008_nu3e-4.log")
runs = (("倍格子 ε 固定（6dz）", "ns_w008_nu3e-4_hires_e6.log"), ("倍格子 ε 半分（3dz）", "ns_w008_nu3e-4_hires.log"), ("種 2", "ns_w008_nu3e-4_s2.log"))
for name, f in runs:
    d = load(f)
    ts = sorted(t for t in d if t in base and t <= 200)
    if not ts:
        continue
    de = max(abs(d[t][0] - base[t][0]) for t in ts)
    dc = max((abs(d[t][1] - base[t][1]) for t in ts if t >= 30), default=float("nan"))
    print(f"{name}: t≤{ts[-1]:.0f}  max|ΔE/E0| {de:.5f}  max|Δ芯 r=0.2|（t≥30） {dc:.4f}")
    for t in (30, 50, 100, 150, 200):
        if t in d and t in base:
            print(f"    t={t}: 芯 {d[t][1]:.3f} 対 {base[t][1]:.3f}  E/E0 {d[t][0]:.5f} 対 {base[t][0]:.5f}")
