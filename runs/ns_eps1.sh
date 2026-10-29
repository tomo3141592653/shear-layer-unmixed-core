cd "$(dirname "$0")/.."
until grep -q ALLDONE ns_w008_nu3e-4_hires.log 2>/dev/null; do sleep 60; done
for nu in 1e-4 2e-4 3e-4 1e-3; do python ns_band.py --w 0.08 --nu $nu --eps 1 --tmax 400 --out ns_w008_nu${nu}_e1 > ns_w008_nu${nu}_e1.log 2>&1; done
echo ALLDONE >> ns_w008_nu1e-3_e1.log
