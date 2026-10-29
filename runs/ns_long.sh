cd "$(dirname "$0")/.."
until grep -q ALLDONE ns_w008_hyper.log 2>/dev/null; do sleep 30; done
python ns_band.py --w 0.08 --nu 1e-4 --tmax 700 --out ns_w008_nu1e-4_long > ns_w008_nu1e-4_long.log 2>&1
echo ALLDONE >> ns_w008_nu1e-4_long.log
