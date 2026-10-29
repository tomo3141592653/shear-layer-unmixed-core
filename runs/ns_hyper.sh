#!/bin/bash
cd "$(dirname "$0")/.."
until grep -q ALLDONE ns_w008_nu1e-3.log 2>/dev/null; do sleep 30; done
python ns_band.py --w 0.08 --nu 0 --nu4 1e-18 --tmax 400 --out ns_w008_hyper > ns_w008_hyper.log 2>&1
echo ALLDONE >> ns_w008_hyper.log
