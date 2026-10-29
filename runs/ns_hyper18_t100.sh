#!/bin/bash
# 本番の超粘性 ν4=1e-18（512×1024、ε=3dz）を t=100 まで回し直し、場を残す（ns_w008_hyper の npz が残っていないため。エンストロフィーと渦度の面積分布を倍の格子と比べる）
cd "$(dirname "$0")/.."
python ns_band.py --w 0.08 --nu 0 --nu4 1e-18 --tmax 100 --out ns_w008_hyper18_t100 > ns_w008_hyper18_t100.log 2>&1
