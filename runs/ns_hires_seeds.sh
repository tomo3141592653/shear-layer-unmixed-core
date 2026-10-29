#!/bin/bash
# 倍の格子（1024×2048、ε=6dz、格子点ごとのゆらぎ）の種 2・3（t≤200）。種1（hires_e6）の 0.473 が分布の裾かを見る
cd "$(dirname "$0")/.."
for s in 2 3; do
  python ns_band.py --w 0.08 --nu 3e-4 --nth 1024 --nz 2048 --dt 0.0025 --eps 6 --seed $s --tmax 200 --out ns_w008_nu3e-4_hires_e6_s$s > ns_w008_nu3e-4_hires_e6_s$s.log 2>&1
done
