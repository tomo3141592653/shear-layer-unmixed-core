"""粘性ありの NS（ns_band.py）の渦度場を、omegapsi2.py と同じ 64×65 格子・同じ壁・同じ循環の規格化で粗視化し、
同じ ψ のビンで ω の条件付き平均をとって、MRS・コアハローの予言との差（RMS）を出す。粘性なし MD（run_w008）と並べる。

    uv run python ns_omegapsi.py 300 400
"""
import math
import os
import sys

import numpy as np


NTH = NZ = 64; ZW = 4.0; w = 0.08
dth = 2 * math.pi / NTH; dz = 2 * ZW / NZ
j = np.arange(1, NZ)
S = np.sin(math.pi * j[:, None] * j[None, :] / NZ) * math.sqrt(2.0 / NZ)
lam_z = (2 * np.cos(math.pi * j / NZ) - 2) / dz**2
k = np.fft.fftfreq(NTH, d=1.0 / NTH)
lam_t = (2 * np.cos(2 * math.pi * k / NTH) - 2) / dth**2
den = lam_t[:, None] + lam_z[None, :]


def solve(G):
    G = G / (G.sum() * dth * dz) * 4 * math.pi
    psi = np.zeros_like(G)
    psi[:, 1:-1] = np.fft.ifft((-(np.fft.fft(G[:, 1:-1], axis=0) @ S) / den) @ S, axis=0).real
    return psi[:, 1:-1].ravel(), G[:, 1:-1].ravel() * w


def coarse_theory(Gf):
    """理論の格子（128×129 でも 256×257 でも）を、半分にするのを繰り返して 64×65 へ。"""
    while Gf.shape[0] > NTH:
        Gt = 0.5 * (Gf[0::2] + Gf[1::2])
        pad = np.concatenate([Gt[:, :1] * 0, Gt, Gt[:, :1] * 0], axis=1)
        Gf = 0.25 * pad[:, 0:-2:2] + 0.5 * pad[:, 1:-1:2] + 0.25 * pad[:, 2::2]
    return Gf


def from_grid(om):
    """ns_band の保存場（θ 256 点 × z 256 点、z=-4..4−Δ）を、64×65 の節点まわりの升目へ平均（循環は solve で規格化）。"""
    nth, nz = om.shape
    th = (np.arange(nth) + 0.0) * 2 * math.pi / nth
    zc = -ZW + np.arange(nz) * 2 * ZW / nz
    TH, ZC = np.meshgrid(th, zc, indexing="ij")
    H, _, _ = np.histogram2d(TH.ravel(), ZC.ravel(), bins=[NTH, NZ + 1],
                             range=[[0, 2 * math.pi], [-ZW - dz / 2, ZW + dz / 2]], weights=om.ravel())
    return H


def from_md(th, zz):
    H, _, _ = np.histogram2d(th % (2 * math.pi), zz, bins=[NTH, NZ + 1], range=[[0, 2 * math.pi], [-ZW - dz / 2, ZW + dz / 2]])
    return H.astype(float)


times = [int(x) for x in sys.argv[1:]] or [300, 400]
fields = {}
md = np.load("run_w008.npz")
fields["粘性なし MD"] = [solve(from_md(*md[f"t{t}"])) for t in times]
for nu in ("1e-4", "3e-4", "1e-3", "hyper"):   # hyper = 超粘性 ν4∇⁸（粘性なしの擬スペクトル）
    try:
        d = np.load("ns_w008_hyper.npz" if nu == "hyper" else f"ns_w008_nu{nu}.npz")
    except FileNotFoundError:
        continue
    fields[f"NS ν={nu}"] = [solve(from_grid(d[f"t{t}"])) for t in times if f"t{t}" in d.files]
mfiles = [("コアハロー", "ch_w008.npz"), ("MRS", "mrs_w008.npz")]
mfiles += [(f"MRS E={e}", f"mrs_w008_e{e}.npz") for e in ("0938", "0.95", "0.90", "0.85", "0.80")]   # エネルギーを減らして解いた MRS
models = {n: solve(coarse_theory(np.load(f)["G"])) for n, f in mfiles if os.path.exists(f)}
allp = np.concatenate([np.concatenate([p for p, _ in v]) for v in fields.values() if v] + [m[0] for m in models.values()])
edges = np.linspace(allp.min(), allp.max(), 41)


def cond(p, o):
    i = np.clip(np.digitize(p, edges) - 1, 0, 39)
    return np.array([o[i == b].mean() if np.any(i == b) else np.nan for b in range(40)])


th_c = {n: cond(*m) for n, m in models.items()}
print(f"t={times}  行: 場、列: 予言との RMS（全体／ハロー 0.001<ω<0.3）と、芯側（ψ 上位 5 ビン）の ω の平均")
for name, v in fields.items():
    if not v:
        continue
    P = np.concatenate([p for p, _ in v]); O = np.concatenate([o for _, o in v])
    c = cond(P, O); out = []
    for n, tc in th_c.items():
        ok = ~np.isnan(c) & ~np.isnan(tc); halo = ok & (c > 0.001) & (c < 0.3)
        out.append(f"{n} {np.sqrt(np.mean((tc[ok]-c[ok])**2)):.3f}/{np.sqrt(np.mean((tc[halo]-c[halo])**2)):.3f}")
    top = c[~np.isnan(c)][-5:]
    print(f"{name:12s} " + "  ".join(out) + f"  芯側 ω {np.mean(top):.3f}")
for n, tc in th_c.items():
    print(f"予言 {n}: 芯側 ω {np.mean(tc[~np.isnan(tc)][-5:]):.3f}")
