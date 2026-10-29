cd "$(dirname "$0")/.."
until grep -q ALLDONE ns_w008_nu1e-3_e1.log 2>/dev/null; do sleep 60; done
python ns_band.py --w 0.08 --nu 0 --nu4 1e-16 --eps 1 --tmax 400 --out ns_w008_hyper16_e1 > ns_w008_hyper16_e1.log 2>&1
echo ALLDONE >> ns_w008_hyper16_e1.log
python ns_band.py --w 0.08 --nu 3e-4 --eps 1.5 --tmax 400 --out ns_w008_nu3e-4_e15 > ns_w008_nu3e-4_e15.log 2>&1
python vic.py --w 0.04 --zw 2 --nth 512 --nz 256 --tmax 400 --out run_w004_zw2 > run_w004_zw2.log 2>&1
echo ALLDONE >> run_w004_zw2.log
