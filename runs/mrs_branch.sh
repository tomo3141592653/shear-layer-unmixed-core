cd "$(dirname "$0")/.."
until [ -f ns_classic.done ]; do sleep 60; done
for b in core zonal two; do python mrs.py --w 0.08 --init $b --out mrs_w008_br_$b > mrs_w008_br_$b.log 2>&1; done
echo ALLDONE >> mrs_w008_br_two.log
