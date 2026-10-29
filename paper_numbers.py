"""論文（paper_draft.md・paper_en.md）の主な数字を、元の計算から丸めずに出し直す。投稿・改稿の前に通して本文と見比べる。
本文の数字をメモの一部の時刻や丸めた表示から書くと、範囲の端がずれる（2026-10-29 に 4 件見つかった）。

    python paper_numbers.py
"""
import re, subprocess, sys
import numpy as np

PY = sys.executable
RX = r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)"


def run(*args):
    return subprocess.run([PY, *args], capture_output=True, text=True).stdout


def nslog(f, tmin=0):
    r = []
    for l in open(f):
        m = re.search(RX, l)
        if m and float(m.group(1)) >= tmin:
            r.append((float(m.group(1)), float(m.group(2)), float(m.group(5))))
    return r


def ratio_from(lam, x):
    p = lam * (lam - 1) / ((lam + 1) * (lam**2 + 1))
    lo, hi = 1.0, 3.0                                  # 式から出る縦横比
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if mid * (mid - 1) / ((mid + 1) * (mid**2 + 1)) < x else (lo, mid)
    return x / p, lo - lam


print("== 3.2 拡散だけの模型との最大差（diffuse_only.py を 6 桁で。NS 側の芯はログが 3 桁なので ±0.0005）")
src = open("diffuse_only.py").read().replace("{E/Eref*e0:.4f} {c:.3f}", "{E/Eref*e0:.6f} {c:.5f}").replace("{ns[0]:.4f} {ns[1]:.3f}", "{ns[0]:.6f} {ns[1]:.5f}")
open("/tmp/_diffuse_hi.py", "w").write(src)
for nu, f, t0, tmax in (("1e-4", "ns_w008_nu1e-4_long", 50, 700), ("3e-4", "ns_w008_nu3e-4", 50, 400), ("1e-4", "ns_w008_nu1e-4_long", 100, 700),
                        ("1e-4", "ns_w008_nu1e-4_e1", 50, 400), ("2e-4", "ns_w008_nu2e-4_e1", 50, 400), ("3e-4", "ns_w008_nu3e-4_e1", 50, 400)):
    out = subprocess.run([PY, "/tmp/_diffuse_hi.py", "--nu", nu, "--file", f + ".npz", "--log", f + ".log", "--t0", str(t0)], capture_output=True, text=True, cwd=".").stdout
    me = mc = 0.0
    for l in out.splitlines():
        if l.startswith("t=") and "NS" in l and int(l.split()[1]) <= tmax:
            x = [float(v) for v in re.findall(r"\d+\.\d+", l)]
            me, mc = max(me, abs(x[0] - x[2])), max(mc, abs(x[1] - x[3]))
    print(f"  {f} t0={t0} t≤{tmax}: |ΔE/E0| ≤ {me:.5f}  |Δ芯| ≤ {mc:.5f}")

print("== 3.2 芯が多段 MRS(E) を横切る点（E で線形に内挿）")
for suf, le in (("", "eps3"), ("_e1", "eps1")):
    g = {"__name__": "x"}
    sys.argv = ["x", suf] if suf else ["x"]
    exec(open("fig_core_energy.py").read().split("fig, ax")[0], g)
    eL = np.array([e for e, _ in g["mrsL"]]); yL = np.array([y for _, y in g["mrsL"]])
    for nu in ("1e-4", "2e-4", "3e-4", "1e-3"):
        f = "ns_w008_nu1e-4_long.log" if (nu == "1e-4" and not suf) else f"ns_w008_nu{nu}{suf}.log"
        d = [(t, E, y, y - np.interp(E, eL, yL)) for t, E, y in nslog(f, 40) if eL.min() <= E <= eL.max()]
        cr = [(a[0] + (b[0] - a[0]) * a[3] / (a[3] - b[3]), a[1] + (b[1] - a[1]) * a[3] / (a[3] - b[3]), a[2])
              for a, b in zip(d, d[1:]) if a[3] > 0 >= b[3]]
        s = " ".join(f"t={t:.0f} E={E:.3f} 芯={y:.3f}" for t, E, y in cr) or "なし"
        print(f"  {le} ν={nu}: 横切り {s}  最初の差 {d[0][3]:+.3f}  最後の差 {d[-1][3]:+.3f}（E={d[-1][1]:.3f}）")

print("== 3.4 表1（classic_table_e.py multi）")
print("\n".join("  " + l for l in run("classic_table_e.py", "multi").splitlines() if "ν=" in l or "MRS 芯" in l))

print("== 3.2 ν→0（超粘性との差/ν の 2 本の比と、2 点の直線を ν=0 へ延ばした値）")
for suf in ("", "_e1"):
    for t in (30, 50, 100, 200):
        def core(f):
            return next(y for tt, _, y in nslog(f) if tt == t)
        y0, y1, y2 = core(f"ns_w008_hyper16{suf}.log"), core(f"ns_w008_nu1e-4{suf}.log"), core(f"ns_w008_nu2e-4{suf}.log")
        print(f"  {suf or 'eps3'} t={t}: 比 {((y0 - y2) / 2) / (y0 - y1):.3f}  延長 {2 * y1 - y2:.3f} 対 超粘性 {y0:.3f}（差 {abs(2 * y1 - y2 - y0):.3f}）")

print("== 4 楕円の芯: 一様パッチの式に対する e_tot/ω_c の比と、式から出る λ との差（t≥100）")
for q in ("0.5", "0.7", "0.3"):
    for f in ("ns_w008_nu1e-4_long", "ns_w008_nu3e-4"):
        rs, ds, lam, t = [], [], None, None
        for l in run("core_strain.py", "--file", f + ".npz", "--q", q).splitlines():
            m = re.match(r"t(\d+):.*縦横比 λ=([\d.]+)", l)
            if m:
                t, lam = int(m.group(1)), float(m.group(2))
            m = re.search(r"e_tot/ω_c=([\d.]+)", l)
            if m and lam and t >= 100:
                r, dl = ratio_from(lam, float(m.group(1))); rs.append(r); ds.append(dl)
        print(f"  q={q} {f}: 比 {min(rs):.3f}〜{max(rs):.3f}  λ の差 {min(ds):.3f}〜{max(ds):.3f}")

print("== 3.1 粘性なし VIC の芯（r=0.2）と、多段の合体（w=0.04、r=0.1、合体後 t≥18）")
v = [float(m.group(1)) for m in re.finditer(r"core\(r=.05/.1/.2\)=[\d.]+/[\d.]+/([\d.]+)", run("vic_core_circle.py"))]
print(f"  VIC w=0.08: {min(v):.3f}〜{max(v):.3f}")
w = [float(m.group(2)) for m in re.finditer(r"t=\s*(\d+) .*= [\d.]+/[\d.]+/([\d.]+)", run("merge_w004.py", "run_w004_zw2.npz", "0.04"))
     if int(m.group(1)) >= 18]
print(f"  VIC w=0.04 合体後: {min(w):.3f}〜{max(w):.3f}")

print("== 2.6 種を変えた差（ν=3e-4、種 1 と 2、t≤400）")
a = {t: (E, y) for t, E, y in nslog("ns_w008_nu3e-4.log")}
b = {t: (E, y) for t, E, y in nslog("ns_w008_nu3e-4_s2.log")}
ts = [t for t in a if t in b]
print(f"  芯（t≥30） {max(abs(a[t][1] - b[t][1]) for t in ts if t >= 30):.4f}  E/E0 {max(abs(a[t][0] - b[t][0]) for t in ts):.5f}")

print("== 3.1 VIC の種 3・4・5（t=35〜100、r=0.2）")
v = [float(m.group(1)) for s in (3, 4, 5) for m in re.finditer(r"core\(r=.05/.1/.2\)=[\d.]+/[\d.]+/([\d.]+)", run("vic_core_circle.py", f"seeds_w008_s{s}.npz"))]
print(f"  {min(v):.3f}〜{max(v):.3f}")

print("== 2.6 縁の物理の幅をそろえた 3 つの格子（ν=3e-4）: 同じ E での芯、多段 MRS(E) との交点、拡散模型との差（t0=50）")
for name, f in (("半分 256×512 ε=1.5dz", "ns_w008_nu3e-4_lores_e15"), ("元 512×1024 ε=3dz", "ns_w008_nu3e-4"), ("倍 1024×2048 ε=6dz", "ns_w008_nu3e-4_hires_e6")):
    r = nslog(f + ".log", 40)
    E = np.array([e for _, e, _ in r]); y = np.array([c for _, _, c in r])
    at = " ".join(f"E={e}:{np.interp(e, E[::-1], y[::-1], left=np.nan, right=np.nan):.3f}" for e in (0.96, 0.95, 0.94, 0.93))
    cr = run("crossing.py", f + ".log").strip().splitlines()[-1].split("横切り ")[1].split("  t=50")[0]
    out = subprocess.run([PY, "/tmp/_diffuse_hi.py", "--nu", "3e-4", "--file", f + ".npz", "--log", f + ".log", "--t0", "50"], capture_output=True, text=True).stdout
    me = mc = 0.0
    for l in out.splitlines():
        if l.startswith("t=") and "NS" in l:
            x = [float(v) for v in re.findall(r"\d+\.\d+", l)]
            me, mc = max(me, abs(x[0] - x[2])), max(mc, abs(x[1] - x[3]))
    print(f"  {name}: {at}  交点 {cr}  拡散模型 |Δ芯|≤{mc:.4f} |ΔE|≤{me:.5f}")
