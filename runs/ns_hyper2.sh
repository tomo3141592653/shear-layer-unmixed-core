cd "$(dirname "$0")/.."
until grep -q ALLDONE mrs_w008_br_two.log 2>/dev/null; do sleep 60; done
python ns_band.py --w 0.08 --nu 0 --nu4 1e-16 --tmax 400 --out ns_w008_hyper16 > ns_w008_hyper16.log 2>&1
echo ALLDONE >> ns_w008_hyper16.log
