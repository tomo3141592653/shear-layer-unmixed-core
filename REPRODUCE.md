---
name: shear-layer-unmixed-core-reproduce
description: Reproduce the figures, Table 1 and the quoted numbers of "An unmixed core in a uniform-vorticity shear layer" from the published result files, or rerun the simulations.
allowed-tools: Bash(git *), Bash(curl *), Bash(tar *), Bash(python *)
---

# Reproduce the shear-layer paper

## 1. Figures, table and numbers from the result files (minutes, CPU only)

```bash
git clone https://github.com/tomo3141592653/shear-layer-unmixed-core
cd shear-layer-unmixed-core
curl -L -o results.tar.gz https://github.com/tomo3141592653/shear-layer-unmixed-core/releases/download/v1.0/shear-layer-results.tar.gz
tar xzf results.tar.gz
pip install numpy matplotlib     # the analysis scripts need only these; torch is for the simulations
python paper_numbers.py            # every number quoted in the text
python fig_core_circle.py          # Figure 1   (add _e1 for the ε = 1dz version)
python fig_core_energy.py          # Figure 2
python fig_core_strain.py          # Figure 3
python fig_nu_limit.py             # Figure 4
python classic_table_e.py multi    # Table 1
```

The figures should match `figures/*.png` exactly.

## 2. Rerun a simulation (GPU, hours)

Each file in `runs/` is the exact command line used, e.g.

```bash
bash runs/ns_sweep.sh              # Navier–Stokes, ν = 1e-4, 3e-4, 1e-3, t ≤ 400
bash runs/vic_main.sh              # inviscid VIC to t = 900 (about 2.3 h on an RTX 2060)
bash runs/ns_conv_fixed_eps.sh     # resolution test at fixed physical edge width
```

`runs/README.md` maps each run to the section of the paper that uses it.

## What to check

- Inviscid core (VIC, r = 0.2): 0.98–1.03γ0 for t = 35–900; seeds 3–5 give 0.94–1.03 for t ≤ 100.
- Viscous core crosses the multi-level MRS(E) once from above (E/E0 ≈ 0.92–0.96 for ν = 3e-4) and moves below it.
- Diffusion-only model from t = 50 reproduces the core within 0.004γ0 and E/E0 within 0.0013 on all grids.
