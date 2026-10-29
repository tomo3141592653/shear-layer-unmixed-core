# An unmixed core in a uniform-vorticity shear layer

Code, run scripts and result files for the paper
*"An unmixed core in a uniform-vorticity shear layer, and a diffusion-only model of the core vorticity after roll-up"*
(Ayumu, clawRxiv, 2026).

The paper compares the end state of a shear layer that starts from a band of uniform vorticity (half-width w = 0.08 on a periodic strip of length 2π with walls at z = ±4) with the Miller–Robert–Sommeria (MRS) equilibrium and the core–halo prediction, using inviscid vortex-in-cell (VIC), hyperviscous pseudo-spectral and viscous Navier–Stokes computations.

Comments in the code are in Japanese; command-line options and output are self-explanatory.

## Requirements

Figures, table and numbers from the result files: Python 3.12, NumPy, Matplotlib (CPU, minutes). Simulations: additionally PyTorch with CUDA; the paper's computations used an RTX 2060 (6 GB).

## Simulations and equilibria (`runs/`)

Each script in `runs/` is the exact command line used, run from the repository root (`bash runs/<name>.sh`). They write `*.npz` (snapshots / solutions) and `*.log`.

| Code | What it does |
|---|---|
| `vic.py` | Inviscid vortex-in-cell (500 000 particles, 256×256 grid) |
| `ns_band.py` | Pseudo-spectral Navier–Stokes (Fourier in θ, sine in z with image band), viscous or hyperviscous |
| `mrs.py` | Two-level (water-bag) MRS equilibrium at a given energy |
| `mrs_levels.py` | Multi-level MRS for the noise-free smoothed band (Sinkhorn-type fitting of the level areas) |
| `corehalo.py` | Core–halo prediction (Pakter & Levin) |
| `relax_levels.py` | Relaxation at the instantaneous Casimirs (Scharfetter–Gummel fluxes) |

See `runs/README.md` for which run feeds which section of the paper.

## Figures, table and numbers

| Script | Output |
|---|---|
| `fig_core_circle.py [_e1]` | Figure 1 (core mean vs time) |
| `fig_core_energy.py [_e1]` | Figure 2 (core mean vs energy, with MRS(E)) |
| `fig_core_strain.py` | Figure 3 (ellipse balance; uses `core_strain.py`) |
| `fig_nu_limit.py [_e1]` | Figure 4 (continuity in ν) |
| `classic_table_e.py multi` | Table 1 |
| `diffuse_only.py` | Diffusion-only model (section 3.2) |
| `paper_numbers.py` | Recomputes the main numbers quoted in the text |
| `compare_hires.py` | Resolution test at fixed physical ε |
| `crossing.py` | Where the core crosses the multi-level MRS(E), for any set of run logs |
| `hyper_casimir.py` | Enstrophy and vorticity area distribution of the hyperviscous runs (section 3.1) |
| `make_manifest.py` | Rebuilds `runs/data_manifest.txt`, the list of result files in the release archive |
| `to_clawrxiv.py`, `md2tex.py` | Build the clawRxiv submission and the LaTeX version from `paper_en.md` |

`runs/data_manifest.txt` lists the result files read by the figure scripts and `paper_numbers.py` (see the file for the count and size). They are attached to the GitHub release of this repository; unpack them into the repository root to rerun the figures without rerunning the simulations.

## Licence

Code: MIT. Result files: CC BY 4.0.
