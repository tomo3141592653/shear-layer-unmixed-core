cd "$(dirname "$0")/.."
PY=python
$PY relax_levels.py --snap ns_w008_nu3e-4.npz --t 300 --K 10 --steps 40000 --sg --every 2000 --out rlx_sg_t300 > rlx_sg_t300.log 2>&1
$PY relax_levels.py --vic run_w008.npz --t 900 --steps 40000 --sg --every 4000 --out rlx_sg_vic > rlx_sg_vic.log 2>&1
for t in 100 200 400; do
$PY relax_levels.py --snap ns_w008_nu3e-4.npz --t $t --K 10 --steps 40000 --sg --every 2000 --out rlx_sg_t$t > rlx_sg_t$t.log 2>&1
done
echo ALLDONE >> rlx_sg_t400.log
for K in 5 20; do
$PY relax_levels.py --snap ns_w008_nu3e-4.npz --t 300 --K $K --steps 40000 --sg --every 2000 --out rlx_sg_t300_K$K > rlx_sg_t300_K$K.log 2>&1
done
