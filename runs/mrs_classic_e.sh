cd "$(dirname "$0")/.."
for w in 0.107 0.214; do for e in 0.95 0.90 0.85; do
python mrs.py --w $w --efrac $e --out mrs_w${w}_e$e > mrs_w${w}_e$e.log 2>&1
done; done
echo ALLDONE >> mrs_w0.214_e0.85.log
