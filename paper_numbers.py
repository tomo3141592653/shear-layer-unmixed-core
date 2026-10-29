"""論文（paper_draft.md・paper_en.md）の主な数字を、元の計算から丸めずに出し直す。投稿・改稿の前に通して本文と見比べる。
本文の数字をメモの一部の時刻や丸めた表示から書くと、範囲の端がずれる（2026-10-29 に 4 件見つかった）。

    python paper_numbers.py
"""
import os, re, subprocess, sys
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
grids = (("半分 256×512 ε=1.5dz", "ns_w008_nu3e-4_lores_e15"), ("元 512×1024 ε=3dz", "ns_w008_nu3e-4"), ("倍 1024×2048 ε=6dz", "ns_w008_nu3e-4_hires_e6"),
         # ゆらぎを 256×512 で作って補間した初期場（--noise-grid 256,512）。半分の格子は lores_e15 と同じ初期場になる
         ("ゆらぎ共通 元 512×1024", "ns_w008_nu3e-4_pn512"), ("ゆらぎ共通 倍 1024×2048", "ns_w008_nu3e-4_pn1024"))
for name, f in grids:
    if not os.path.exists(f + ".npz"):
        print(f"  {name}: {f}.npz が無い")
        continue
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

def core_at(f, Es=(0.96, 0.95, 0.94, 0.93), tmax=200):
    r = [x for x in nslog(f + ".log", 40) if x[0] <= tmax]
    E = np.array([e for _, e, _ in r]); y = np.array([c for _, _, c in r])
    return [float(np.interp(e, E[::-1], y[::-1], left=np.nan, right=np.nan)) for e in Es]

print("== 2.6 同じ初期場の 3 格子の差（t≥40、芯 r=0.2 と E/E0）")
logs = {k: {t: (E, y) for t, E, y in nslog(f"ns_w008_nu3e-4_{k}.log", 40)} for k in ("lores_e15", "pn512", "pn1024") if os.path.exists(f"ns_w008_nu3e-4_{k}.log")}
for a_, b_ in (("lores_e15", "pn512"), ("pn512", "pn1024")):
    if a_ in logs and b_ in logs:
        ts = [t for t in logs[a_] if t in logs[b_]]
        print(f"  {a_}↔{b_}: t={min(ts):.0f}〜{max(ts):.0f} 芯 {max(abs(logs[a_][t][1] - logs[b_][t][1]) for t in ts):.3f}  E/E0 {max(abs(logs[a_][t][0] - logs[b_][t][0]) for t in ts):.5f}")

print("== 2.6 半分の格子の種 1〜6（E/E0=0.96/0.95/0.94/0.93 の芯、交点）")
rows, crs = [], []
for s in ("", "_s2", "_s3", "_s4", "_s5", "_s6"):
    f = "ns_w008_nu3e-4_lores_e15" + s
    rows.append(core_at(f))
    crs.append(float(re.search(r"横切り t=\d+ E=([\d.]+)", run("crossing.py", f + ".log")).group(1)))
a = np.array(rows)
print("  平均 " + " ".join(f"{x:.3f}" for x in a.mean(0)) + "  標準偏差 " + " ".join(f"{x:.3f}" for x in a.std(0, ddof=1)) + "  幅 " + " ".join(f"{x:.3f}〜{y:.3f}" for x, y in zip(a.min(0), a.max(0))))
print(f"  交点 E/E0={min(crs):.3f}〜{max(crs):.3f}")
print("== 2.6 ゆらぎの振幅だけを変える（半分の格子・種1、E/E0=0.95 の芯）")
base = core_at("ns_w008_nu3e-4_lores_e15")[1]
for A in ("5e-4", "2.5e-4"):
    print(f"  振幅 {A}: {core_at('ns_w008_nu3e-4_lores_a' + A)[1]:.3f}（元 {base:.3f}、差 {core_at('ns_w008_nu3e-4_lores_a' + A)[1] - base:+.3f}）")
for s in ("_s2", "_s3"):
    if os.path.exists(f"ns_w008_nu3e-4_lores_a2.5e-4{s}.log"):
        v, b0 = core_at("ns_w008_nu3e-4_lores_a2.5e-4" + s)[1], core_at("ns_w008_nu3e-4_lores_e15" + s)[1]
        print(f"  振幅 2.5e-4 種{s[2:]}: {v:.3f}（元 {b0:.3f}、差 {v - b0:+.3f}）")
print("== 2.6 格子点ごとのゆらぎの計算の E/E0=0.95 の芯（半分の格子の種 1〜6 の平均との差）")
m = a.mean(0)[1]
for name, f in (("元 種1", "ns_w008_nu3e-4"), ("元 種2", "ns_w008_nu3e-4_s2"), ("倍 種1", "ns_w008_nu3e-4_hires_e6"), ("元 細かい成分を落とした", "ns_w008_nu3e-4_lp256")):
    if os.path.exists(f + ".log"):
        v = core_at(f, tmax=400)[1]
        print(f"  {name}: {v:.3f}（平均 {m:.3f} との差 {v - m:+.3f}）")
print("== 2.6 元の格子・格子点ごとのゆらぎの種 1〜5（E/E0=0.96/0.95/0.94/0.93 の芯、交点、拡散模型）")
rows, crs = [], []
for s in ("", "_s2", "_s3", "_s4", "_s5"):
    f = "ns_w008_nu3e-4" + s
    if not os.path.exists(f + ".log"):
        continue
    rows.append(core_at(f))
    crs.append(float(re.search(r"横切り t=\d+ E=([\d.]+)", run("crossing.py", f + ".log")).group(1)))
b = np.array(rows)
print(f"  {len(rows)} 種 平均 " + " ".join(f"{x:.3f}" for x in b.mean(0)) + "  標準偏差 " + " ".join(f"{x:.3f}" for x in b.std(0, ddof=1)) + f"  半分の格子の平均との差（E=0.95） {b.mean(0)[1] - a.mean(0)[1]:+.3f}")
print(f"  交点 E/E0={min(crs):.3f}〜{max(crs):.3f}")
print("== 2.6 倍の格子・格子点ごとのゆらぎの種 1〜3（E/E0=0.95 の芯、交点、拡散模型 t0=50）")
for s in ("", "_s2", "_s3"):
    f = "ns_w008_nu3e-4_hires_e6" + s
    if not os.path.exists(f + ".npz"):
        continue
    cr = re.search(r"横切り t=\d+ E=([\d.]+)", run("crossing.py", f + ".log")).group(1)
    out = subprocess.run([PY, "/tmp/_diffuse_hi.py", "--nu", "3e-4", "--file", f + ".npz", "--log", f + ".log", "--t0", "50"], capture_output=True, text=True).stdout
    me = mc = 0.0
    for l in out.splitlines():
        if l.startswith("t=") and "NS" in l:
            x = [float(v) for v in re.findall(r"\d+\.\d+", l)]
            me, mc = max(me, abs(x[0] - x[2])), max(mc, abs(x[1] - x[3]))
    print(f"  種{s[2:] or '1'}: {core_at(f)[1]:.3f}  交点 E/E0={cr}  拡散模型 |Δ芯|≤{mc:.4f} |ΔE|≤{me:.5f}")
print("== 2.6 倍の格子・種1 から元の格子で表せない成分を落とす（lp512）、E/E0=0.97/0.96/0.95 の芯")
if os.path.exists("ns_w008_nu3e-4_hires_lp512.log"):
    for f in ("ns_w008_nu3e-4_hires_e6", "ns_w008_nu3e-4_hires_lp512"):
        print(f"  {f}: " + " ".join(f"{v:.3f}" for v in core_at(f, Es=(0.97, 0.96, 0.95))))
from math import sqrt
def welch(x, y):
    x, y = np.array(x), np.array(y)
    return (y.mean() - x.mean()) / sqrt(x.var(ddof=1) / len(x) + y.var(ddof=1) / len(y))
hs = [core_at("ns_w008_nu3e-4_hires_e6" + s)[1] for s in ("", "_s2", "_s3") if os.path.exists("ns_w008_nu3e-4_hires_e6" + s + ".log")]
print(f"  倍の格子 {len(hs)} 種 平均 {np.mean(hs):.3f} 標準偏差 {np.std(hs, ddof=1):.3f}  ウェルチ t: 元−半分 {welch(a[:, 1], b[:, 1]):.2f}  倍−元 {welch(b[:, 1], hs):.2f}")
print("== 3.1 VIC の収束（粒子 10⁶・格子 256/384、t=35/48/100 の芯 r=0.05/0.1/0.2）")
for f in ("run_w008.npz", "vic_w008_np1m.npz", "vic_w008_np1m_g384.npz"):
    if os.path.exists(f):
        v = [(m.group(1), m.group(2), m.group(3), m.group(4)) for m in re.finditer(r"t=(\d+) core\(r=.05/.1/.2\)=([\d.]+)/([\d.]+)/([\d.]+)", run("vic_core_circle.py", f)) if m.group(1) in ("35", "48", "100")]
        print(f"  {f}: " + "  ".join(f"t={t} {a}/{b}/{c}" for t, a, b, c in v))
print("== 3.1 超粘性の収束（倍の格子 1024×2048・ν4=3.9e-21 と 512×1024・ν4=1e-18、t=30〜100 の芯 r=0.2 と中心 r=0.05）")
if os.path.exists("ns_w008_hyper_g1024.log"):
    def hl(f):
        return {float(t): (float(a), float(c), float(e)) for t, e, a, b, c in re.findall(r"t=\s*([\d.]+) E/E0=([\d.]+) core\(r=\.05/\.1/\.2\)=([\d.]+)/([\d.]+)/([\d.]+)", open(f).read())}
    A, G = hl("ns_w008_hyper18_t100.log"), hl("ns_w008_hyper_g1024.log")   # 元の格子は t=100 まで回し直した版（ns_w008_hyper.log と時間の欄以外同じ）
    ts = [t for t in G if 30 <= t <= 100]
    print(f"  r=0.2: 倍の格子 {min(G[t][1] for t in ts):.3f}〜{max(G[t][1] for t in ts):.3f}  最大差 {max(abs(G[t][1] - A[t][1]) for t in ts):.3f}  中心 t=100 {G[100.0][0]:.3f}  |E/E0−1|≤{max(abs(G[t][2] - 1) for t in G):.5f}")
print("== 3.1 超粘性のエンストロフィーと渦度の面積分布（hyper_casimir.py、元 ν4=1e-18 と倍 ν4=3.9e-21）、ν4 k_c^8")
print(f"  ν4 k_c^8: 元 {1e-18 * (512 / 3)**8:.2f}  倍 {3.9e-21 * (1024 / 3)**8:.2f}")
if os.path.exists("ns_w008_hyper18_t100.npz"):
    print("  " + run("hyper_casimir.py", "ns_w008_hyper18_t100.npz", "ns_w008_hyper_g1024.npz").replace("\n", "\n  ").rstrip())
print("== 4 議論: 芯の中心が MRS の最大渦度の何倍か（VIC の中心 1.0、超粘性の中心 0.96 / γ0）")
m2, mL = (float(np.load(f)["G"].max()) * 0.08 for f in ("mrs_w008_256.npz", "mrsL_eps3.npz"))
print(f"  2値 MRS の最大 {m2:.3f} → VIC {1.0 / m2:.2f} 倍   多段 MRS の最大 {mL:.3f} → 超粘性 {0.96 / mL:.2f} 倍")
print("== 3.3 その時刻の Casimir を保つ緩和（relax_levels.py --sg、最初と最後の S と芯 r=0.2）")
for f in ("rlx_sg_t100", "rlx_sg_t200", "rlx_sg_t300", "rlx_sg_t400", "rlx_sg_t300_K5", "rlx_sg_t300_K20", "rlx_sg_vic", "rlx_sg_mrsstart"):
    if os.path.exists(f + ".log"):
        st = re.findall(r"step=(\d+) .*?S=([-\d.]+) .*?芯 r=0.2 ([\d.]+)", open(f + ".log").read())
        if st:
            print(f"  {f}: S {float(st[0][1]):.2f}→{float(st[-1][1]):.2f}  芯 {st[0][2]}→{st[-1][2]}（{st[-1][0]} 歩）")
print("== 2.4 初期の帯で中間の渦度（0.05〜0.95γ0）が運ぶ循環の割合（連続の tanh 分布／NS の格子点）")
for k in (3, 1):
    eps = k * 16 / 1024
    for nm, zz in (("連続", np.linspace(-4, 4, 800001)), ("格子点", -4 + np.arange(512) * 16 / 1024)):
        s = 0.5 * (np.tanh((0.08 - np.abs(zz)) / eps) + 1)
        print(f"  ε={k}dz {nm}: {s[(s > 0.05) & (s < 0.95)].sum() / s.sum() * 100:.1f}%", end="")
    print()
