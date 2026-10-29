#!/bin/bash
# ε を物理的に固定した収束試験: 倍の格子（1024×2048）で縁を ε=6dz にし、元の格子（512×1024）の ε=3dz と同じ物理の幅にする。
# ns_conv.sh の倍格子は ε=3dz（物理の幅が半分）だったので、離散化誤差と ε の変化が分けられなかった（査読3回目の指摘）。
cd "$(dirname "$0")/.."
python ns_band.py --w 0.08 --nu 3e-4 --nth 1024 --nz 2048 --dt 0.0025 --eps 6 --tmax 200 --out ns_w008_nu3e-4_hires_e6 > ns_w008_nu3e-4_hires_e6.log 2>&1
echo ALLDONE >> ns_w008_nu3e-4_hires_e6.log
# 半分の格子（256×512、ε=1.5dz、dt=0.01）。3 つの格子で収束の速さと行き先を見積もる
python ns_band.py --w 0.08 --nu 3e-4 --nth 256 --nz 512 --dt 0.01 --eps 1.5 --tmax 400 --out ns_w008_nu3e-4_lores_e15 > ns_w008_nu3e-4_lores_e15.log 2>&1
