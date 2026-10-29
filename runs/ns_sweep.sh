#!/bin/bash
cd "$(dirname "$0")/.."
for nu in 1e-4 3e-4 1e-3; do
  python ns_band.py --w 0.08 --nu $nu --tmax 400 --out ns_w008_nu$nu > ns_w008_nu$nu.log 2>&1
done
echo ALLDONE >> ns_w008_nu1e-3.log
