#!/bin/bash
# 初期のゆらぎを物理の尺度でそろえた格子収束（ゆらぎは 256×512 で作ってフーリエで補間、縁の物理の幅も同じ）。
# 256×512 はこの定義で ns_w008_nu3e-4_lores_e15 と完全に同じ初期場になるので回し直さない
cd "$(dirname "$0")/.."
PY=python
$PY ns_band.py --w 0.08 --nu 3e-4 --eps 3 --noise-grid 256,512 --tmax 200 --out ns_w008_nu3e-4_pn512 > ns_w008_nu3e-4_pn512.log 2>&1
$PY ns_band.py --w 0.08 --nu 3e-4 --nth 1024 --nz 2048 --dt 0.0025 --eps 6 --noise-grid 256,512 --tmax 200 --out ns_w008_nu3e-4_pn1024 > ns_w008_nu3e-4_pn1024.log 2>&1
echo ALLDONE > ns_conv_physnoise.done
