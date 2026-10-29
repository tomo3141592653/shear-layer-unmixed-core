#!/bin/bash
# NS ν=2e-4（10/28 21:25 に事前予言を書いてから回した分）。元の起動コマンドは残っていなかったが、
# ns_w008_nu2e-4.log の先頭（grid 512x1024、dt 0.005、ε・ゆらぎは既定）が ns_band.py の既定値と一致する。
cd "$(dirname "$0")/.."
python ns_band.py --w 0.08 --nu 2e-4 --tmax 400 --out ns_w008_nu2e-4 > ns_w008_nu2e-4.log 2>&1
