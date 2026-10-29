#!/bin/bash
# 半分の格子（256×512、ε=1.5dz）で種 2・3。格子による芯の差（最大 0.05γ0、単調でない）が、ゆらぎの実現の違いでどこまで説明できるかを見る
cd "$(dirname "$0")/.."
for s in 2 3; do
  python ns_band.py --w 0.08 --nu 3e-4 --nth 256 --nz 512 --dt 0.01 --eps 1.5 --seed $s --tmax 400 --out ns_w008_nu3e-4_lores_e15_s$s > ns_w008_nu3e-4_lores_e15_s$s.log 2>&1
done
echo ALLDONE > ns_lores_seeds.done
