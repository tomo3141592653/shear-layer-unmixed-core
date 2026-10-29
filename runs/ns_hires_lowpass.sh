#!/bin/bash
# 倍の格子（1024×2048、ε=6dz）・種1 の格子点ごとのゆらぎから、元の格子（512×1024）で表せない成分だけを落とす
cd "$(dirname "$0")/.."
python ns_band.py --w 0.08 --nu 3e-4 --nth 1024 --nz 2048 --dt 0.0025 --eps 6 --noise-lowpass 512,1024 --tmax 200 --out ns_w008_nu3e-4_hires_lp512 > ns_w008_nu3e-4_hires_lp512.log 2>&1
