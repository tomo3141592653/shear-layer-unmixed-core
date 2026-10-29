#!/bin/bash
# 半分の格子の種 4〜6（t≤200）。乱数の模様による芯の値の幅を、種 6 本で見る
cd "$(dirname "$0")/.."
for s in 4 5 6; do
  python ns_band.py --w 0.08 --nu 3e-4 --nth 256 --nz 512 --dt 0.01 --eps 1.5 --seed $s --tmax 200 --out ns_w008_nu3e-4_lores_e15_s$s > ns_w008_nu3e-4_lores_e15_s$s.log 2>&1
done
