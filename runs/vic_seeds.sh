#!/bin/bash
# 粘性なし VIC の種 3・4・5（w=0.08、t=100 まで、t=60 からは 2 刻みでも保存）。10/19 に別の問い（合体と島）のために回したもの。
# 起動コマンドは残っていなかったが、notes.md 10/19 01:34 と、同じ夜の w=0.13 の種（vic_track.py --tmax 100 --fine 2 --fine-from 60）と同じ形。
cd "$(dirname "$0")/.."
for s in 3 4 5; do
  python vic_track.py --w 0.08 --seed $s --tmax 100 --fine 2 --fine-from 60 --out seeds_w008_s$s > seeds_w008_s$s.log 2>&1
done
