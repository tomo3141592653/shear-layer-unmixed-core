# 計算の起動コマンド

論文（paper_draft.md・paper_en.md）の計算をどう回したか。どれも `bash runs/<名前>.sh` で、研究フォルダに結果（.npz）とログ（.log）を書く。
もとはセッションの一時フォルダに置いていたもので、2026-10-29 にここへ移した。

| 列 | 中身 | 論文で使うところ |
|---|---|---|
| theory_main.sh | コアハローと 2 値 MRS（E0、128 格子）、2 値 MRS（256 格子、E0 と 0.938） | 2.4、3.1、図1・2 |
| vic_main.sh | 粘性なし VIC、w=0.08、t=900 | 3.1、図1、表1 |
| vic_seeds.sh | 粘性なし VIC の種 3・4・5、t=100 まで（vic_track.py） | 3.1 |
| ns_conv_fixed_eps.sh | 倍の格子で ε=6dz（元の格子の ε=3dz と同じ物理の幅） | 2.6 |
| ns_lores_seeds.sh / ns_lores_seeds2.sh | 半分の格子（256×512、ε=1.5dz）の種 2・3（t=400）と 4〜6（t=200） | 2.6、図2の帯、限界(5) |
| ns_conv_physnoise.sh | ゆらぎを 256×512 で作って補間した同じ初期場で、元の格子と倍の格子（--noise-grid） | 2.6 |
| ns_noise_amp_seeds.sh | 半分の格子で振幅 1/4 を種 2・3 でも | 2.6 |
| ns_noise_amp.sh | 半分の格子・種1で、ゆらぎの振幅だけを 1/2・1/4 に | 2.6 |
| ns_base_seeds.sh | 元の格子・格子点ごとのゆらぎの種 3〜5（t=200） | 2.6、限界(5) |
| ns_hires_seeds.sh | 倍の格子・格子点ごとのゆらぎの種 2・3（t=200） | 2.6、限界(5) |
| ns_hires_lowpass.sh | 倍の格子・種1 のゆらぎから元の格子で表せない成分を落とす | 2.6、限界(5) |
| inviscid_conv.sh | 粘性なしの収束: VIC 粒子 10⁶・格子 256/384、超粘性 1024×2048（t=100） | 3.1、2.6 |
| ns_noise_lowpass.sh | 元の格子・種1の格子点ごとのゆらぎから、半分の格子で表せない成分を落とす（--noise-lowpass） | 2.6 |
| ns_sweep.sh / ns_nu2e-4.sh | 擬スペクトル NS、ν=1e-4・3e-4・1e-3 と 2e-4、ε=3dz、t=400 | 3.2、図1・2・4 |
| ns_long.sh | NS ν=1e-4 を t=700 まで | 3.2、図1・2・3、拡散模型 |
| ns_hyper.sh / ns_hyper2.sh | 超粘性（ν4=1e-18 / 1e-16） | 3.1、図1・4 |
| ns_eps1.sh / ns_hyper_e1.sh | 縁を ε=1dz にした NS と超粘性、ε=1.5dz、多段の合体の VIC（w=0.04） | 2.6、3.2、限界(3) |
| ns_conv.sh | 種 2、倍の格子（1024×2048） | 2.6・2.7 |
| ns_classic.sh | 古典の厚さ（w/L=0.017・0.034）の NS・VIC・2 値 MRS | 3.4、表1 |
| mrs_esweep.sh / mrs_branch.sh / mrs_classic_e.sh | 2 値 MRS のエネルギーを変えた解・枝の確認・古典の厚さ | 3.2・3.4、図2 |
| mrsL.sh / mrsL_classic.sh | 多段 MRS（K=16、512 格子）と E を変えた解 | 2.4、3.2・3.4、図1・2、表1 |
| rlx_sg.sh / rlxK.sh | その時刻の Casimir を保つ緩和（SG 流束、K=5/10/20） | 3.3 |
| rlx_times.sh / rlx_renorm.sh | SG 流束に作り直す前の緩和（格子の分散が芯を混ぜた） | 使わない。経緯の記録 |

図と数字は `fig_*.py` と `classic_table_e.py`、本文の主な数字の照合は `paper_numbers.py`。
