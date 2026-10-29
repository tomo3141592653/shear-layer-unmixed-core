#!/bin/bash
# 半分の格子（256×512, ε=1.5dz）で、ゆらぎの振幅だけを 1/2・1/4 にする。
# 格子点ごとの白色雑音は、点が4倍になると大きなスケールの成分の振幅が半分になる。
# 元の格子・倍の格子の雑音が大きなスケールで持つ振幅をこの格子で再現し、芯の値が振幅でどれだけ動くかを見る
cd "$(dirname "$0")/.."
PY=python
for A in 5e-4 2.5e-4; do
  $PY ns_band.py --w 0.08 --nu 3e-4 --nth 256 --nz 512 --dt 0.01 --eps 1.5 --noise $A --tmax 200 --out ns_w008_nu3e-4_lores_a$A > ns_w008_nu3e-4_lores_a$A.log 2>&1
done
