cd "$(dirname "$0")/.."
for t in 100 200 400; do
python relax_levels.py --snap ns_w008_nu3e-4.npz --t $t --K 10 --steps 30000 --exactE --every 3000 --out rlx_K10_t$t > rlx_K10_t$t.log 2>&1
done
echo ALLDONE >> rlx_K10_t400.log
