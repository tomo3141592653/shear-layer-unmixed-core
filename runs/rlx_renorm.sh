cd "$(dirname "$0")/.."
PY=python
$PY relax_levels.py --vic run_w008.npz --t 900 --steps 40000 --exactE --renorm --every 4000 --out rlx_vic_rn > rlx_vic_rn.log 2>&1
for t in 300 100 200 400; do
$PY relax_levels.py --snap ns_w008_nu3e-4.npz --t $t --K 10 --steps 40000 --exactE --renorm --every 4000 --out rlx_K10rn_t$t > rlx_K10rn_t$t.log 2>&1
done
echo ALLDONE >> rlx_K10rn_t400.log
