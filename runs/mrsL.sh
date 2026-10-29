cd "$(dirname "$0")/.."
PY=python
for e in 3 1; do $PY mrs_levels.py --eps $e --K 16 --nz 512 --out mrsL_eps${e} > mrsL_eps${e}.log 2>&1; done
$PY mrs_levels.py --eps 0.05 --K 2 --nz 512 --out mrsL_check512 > mrsL_check512.log 2>&1
echo ALLDONE >> mrsL_check512.log
for f in 0.95 0.90 0.85; do $PY mrs_levels.py --eps 3 --K 16 --nz 512 --efrac $f --out mrsL_eps3_e$f > mrsL_eps3_e$f.log 2>&1; done
echo ALLDONE >> mrsL_eps3_e0.85.log
for f in 0.95 0.90 0.85; do $PY mrs_levels.py --eps 1 --K 16 --nz 512 --efrac $f --out mrsL_eps1_e$f > mrsL_eps1_e$f.log 2>&1; done
echo ALLDONE >> mrsL_eps1_e0.85.log
