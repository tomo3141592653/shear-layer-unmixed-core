# 古典（Sommeria, Staquet & Robert 1991）の帯の厚さ 0.017L・0.034L（L=2π）で、粘性あり・なし・MRS をそろえる
cd "$(dirname "$0")/.."
PY=python
until grep -q ALLDONE ns_w008_nu1e-4_long.log 2>/dev/null; do sleep 30; done
for w in 0.107 0.214; do
  for nu in 3.3e-4 6.7e-4 1.3e-3; do
    $PY ns_band.py --w $w --nu $nu --tmax 120 --out ns_w${w}_nu$nu > ns_w${w}_nu$nu.log 2>&1
  done
  $PY vic.py --w $w --tmax 120 --out run_w$w > run_w$w.log 2>&1
done
until grep -q ALLDONE mrs_w008_e0.80.log 2>/dev/null; do sleep 30; done
for w in 0.107 0.214; do $PY mrs.py --w $w --out mrs_w$w > mrs_w$w.log 2>&1; done
echo ALLDONE > ns_classic.done
