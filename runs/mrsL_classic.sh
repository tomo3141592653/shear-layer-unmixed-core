cd "$(dirname "$0")/.."
PY=python
while kill -0 1768586 2>/dev/null; do sleep 60; done
for w in 0.107 0.214; do for f in 1.0 0.95 0.90 0.85; do
  $PY mrs_levels.py --w $w --eps 3 --K 16 --nz 512 --efrac $f --out mrsL_w${w}_e$f > mrsL_w${w}_e$f.log 2>&1
done; done
echo ALLDONE >> mrsL_w0.214_e0.85.log
