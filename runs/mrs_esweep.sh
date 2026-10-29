cd "$(dirname "$0")/.."
while kill -0 1131588 2>/dev/null; do sleep 30; done
for e in 0.95 0.90 0.85 0.80; do python mrs.py --w 0.08 --efrac $e --out mrs_w008_e$e > mrs_w008_e$e.log 2>&1; done
echo ALLDONE >> mrs_w008_e0.80.log
