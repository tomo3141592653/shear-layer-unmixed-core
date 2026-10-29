"""公開する結果ファイルの一覧 runs/data_manifest.txt を作り直す。
今の一覧に載っているファイルに、引数で渡した名前（.log と .npz の組なら拡張子なしでよい）を足し、大きさと合計を測り直す。

    python make_manifest.py ns_w008_nu3e-4_pn512 ns_w008_nu3e-4_s3 ...
"""
import os, sys

P = "runs/data_manifest.txt"
names = [l.split("\t")[0] for l in open(P) if l.strip() and not l.startswith("#")]
for a in sys.argv[1:]:
    for f in ([a] if a.endswith((".log", ".npz")) else [a + ".log", a + ".npz"]):
        if not os.path.exists(f):
            sys.exit(f"{f} が無い")
        if f not in names:
            names.append(f)
missing = [f for f in names if not os.path.exists(f)]
if missing:
    sys.exit("一覧にあるのに無い: " + " ".join(missing))
names.sort()
total = sum(os.path.getsize(f) for f in names)
with open(P, "w") as w:
    w.write(f"# 図と paper_numbers.py が開く結果ファイルと、2.6 節・限界(5)で引く種・格子の計算。{len(names)} 個・{total / 2**20:.0f}M\n")
    for f in names:
        w.write(f"{f}\t{os.path.getsize(f) // 1024} KB\n")
print(f"{len(names)} 個・{total / 2**20:.0f}M")
