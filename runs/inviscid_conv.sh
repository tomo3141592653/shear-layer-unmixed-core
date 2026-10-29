#!/bin/bash
# 粘性なしの側の収束（査読5回目の分析の指摘）。w=0.08、t≤100（巻き上がりから芯ができるまで）
# VIC: 粒子 10⁶ で格子 256 と 384（本番は 5×10⁵・256）
# 超粘性: 倍の格子 1024×2048、ε=6dz（縁の物理の幅は本番の 512×1024・ε=3dz と同じ）、ν4 は打ち切り波数での減衰率をそろえて 1e-18/2⁸
cd "$(dirname "$0")/.."
PY=python
$PY vic.py --w 0.08 --np 1000000 --tmax 100 --out vic_w008_np1m > vic_w008_np1m.log 2>&1
$PY vic.py --w 0.08 --np 1000000 --nth 384 --nz 384 --tmax 100 --out vic_w008_np1m_g384 > vic_w008_np1m_g384.log 2>&1
$PY ns_band.py --w 0.08 --nu 0 --nu4 3.9e-21 --nth 1024 --nz 2048 --dt 0.0025 --eps 6 --tmax 100 --out ns_w008_hyper_g1024 > ns_w008_hyper_g1024.log 2>&1
