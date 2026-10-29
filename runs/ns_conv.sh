# 種と格子の確認: ν=3e-4 を種2で、また格子を倍にして（t=200 まで）
cd "$(dirname "$0")/.."
until grep -q ALLDONE ns_w008_hyper16.log 2>/dev/null; do sleep 60; done
python ns_band.py --w 0.08 --nu 3e-4 --seed 2 --tmax 400 --out ns_w008_nu3e-4_s2 > ns_w008_nu3e-4_s2.log 2>&1
python ns_band.py --w 0.08 --nu 3e-4 --nth 1024 --nz 2048 --dt 0.0025 --tmax 200 --out ns_w008_nu3e-4_hires > ns_w008_nu3e-4_hires.log 2>&1
echo ALLDONE >> ns_w008_nu3e-4_hires.log
