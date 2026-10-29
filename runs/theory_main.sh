#!/bin/bash
# 事前の予言: コアハローと 2 値 MRS（E0、128 格子）、2 値 MRS（E0 と E/E0=0.938、256 格子）。
# 元の起動コマンドは残っていなかった。10/29 に今のコードで 128 格子の 2 本を解き直し、保存してある解と G が完全に一致（差 0）することを確かめた。
cd "$(dirname "$0")/.."
PY=python
$PY corehalo.py --w 0.08 --nth 128 --nz 128 --out ch_w008 > ch_w008.log 2>&1
$PY mrs.py --w 0.08 --nth 128 --nz 128 --out mrs_w008 > mrs_w008.log 2>&1
$PY mrs.py --w 0.08 --out mrs_w008_256 > mrs_w008_256.log 2>&1
$PY mrs.py --w 0.08 --efrac 0.938 --out mrs_w008_e0938 > mrs_w008_e0938.log 2>&1
