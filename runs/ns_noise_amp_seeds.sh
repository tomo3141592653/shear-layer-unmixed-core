#!/bin/bash
# 半分の格子で、ゆらぎの振幅 1/4（倍の格子の格子点ごとのゆらぎが大きなスケールで持つ振幅）を種 2・3 でも
cd "$(dirname "$0")/.."
for s in 2 3; do
  python ns_band.py --w 0.08 --nu 3e-4 --nth 256 --nz 512 --dt 0.01 --eps 1.5 --noise 2.5e-4 --seed $s --tmax 200 --out ns_w008_nu3e-4_lores_a2.5e-4_s$s > ns_w008_nu3e-4_lores_a2.5e-4_s$s.log 2>&1
done
