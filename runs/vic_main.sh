#!/bin/bash
# 粘性なし VIC の主計算（w=0.08、粒子 50 万、256×256、dt=0.01、t=900、種 1）。
# 元の起動コマンドは残っていなかったが、run_w008.log の先頭（N=500000）と終わり（t=900）が vic.py の既定値と一致する。
cd "$(dirname "$0")/.."
python vic.py --w 0.08 --out run_w008 > run_w008.log 2>&1
