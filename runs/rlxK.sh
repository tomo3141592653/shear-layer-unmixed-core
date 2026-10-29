cd "$(dirname "$0")/.."
PY=python
while kill -0 1832324 2>/dev/null; do sleep 60; done
for K in 5 20; do
  $PY relax_levels.py --snap ns_w008_nu3e-4.npz --t 300 --K $K --steps 40000 --sg --every 2000 --out rlx_sg_t300_K$K > rlx_sg_t300_K$K.log 2>&1
done
echo ALLDONE >> rlx_sg_t300_K20.log
