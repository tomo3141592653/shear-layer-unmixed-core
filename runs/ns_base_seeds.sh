#!/bin/bash
# 元の格子（512×1024, ε=3dz、格子点ごとのゆらぎ）の種 3〜5（t≤200）。元の格子での種の分布を、半分の格子の種 1〜6 と比べる
cd "$(dirname "$0")/.."
for s in 3 4 5; do
  python ns_band.py --w 0.08 --nu 3e-4 --eps 3 --seed $s --tmax 200 --out ns_w008_nu3e-4_s$s > ns_w008_nu3e-4_s$s.log 2>&1
done
