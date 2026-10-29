# Run scripts

Each script is the command line used for the paper; run from the repository root with `bash runs/<name>.sh`.

| Script | Contents | Used in |
|---|---|---|
| `theory_main.sh` | Core–halo and two-level MRS at E0 (128 grid), two-level MRS at E0 (256 grid) | §2.4, §3.1, Fig. 1 |
| `vic_main.sh` | Inviscid VIC, w = 0.08, t ≤ 900 (seed 1) | §3.1, Fig. 1, Table 1 |
| `vic_seeds.sh` | Inviscid VIC seeds 3–5 to t = 100 (`vic_track.py`, extra snapshots every 2 from t = 60) | §3.1 |
| `ns_sweep.sh`, `ns_nu2e-4.sh` | Navier–Stokes, ν = 1e-4, 2e-4, 3e-4, 1e-3, ε = 3dz, t ≤ 400 | §3.2, Figs 1, 2, 4 |
| `ns_long.sh` | Navier–Stokes ν = 1e-4 to t = 700 | §3.2, Figs 1–3, diffusion-only model |
| `ns_hyper.sh`, `ns_hyper2.sh` | Hyperviscous runs (ν₄ = 1e-18, 1e-16) | §3.1, Figs 1, 4 |
| `ns_eps1.sh`, `ns_hyper_e1.sh` | Sharper edge ε = 1dz (viscous and hyperviscous), ε = 1.5dz, and the multi-merger VIC (w = 0.04, walls ±2) | §2.6, §3.2, Limitations (3) |
| `ns_conv.sh` | Seed 2; doubled grid with ε = 3dz | §2.6 |
| `ns_conv_fixed_eps.sh` | Doubled grid with ε = 6dz (same physical edge width as the base grid) | §2.6 |
| `ns_classic.sh` | Classical band thicknesses (w/L = 0.017, 0.034): Navier–Stokes, VIC, two-level MRS | §3.4, Table 1 |
| `mrs_esweep.sh`, `mrs_branch.sh`, `mrs_classic_e.sh` | Two-level MRS at reduced energies, branch check, classical thicknesses | §3.2, §3.4, Fig. 2 |
| `mrsL.sh`, `mrsL_classic.sh` | Multi-level MRS (K = 16, 512 grid) at E/E0 = 1.00–0.85 | §2.4, §3.2, §3.4, Figs 1, 2, Table 1 |
| `rlx_sg.sh`, `rlxK.sh` | Relaxation at the instantaneous Casimirs (Scharfetter–Gummel fluxes, K = 5, 10, 20) | §3.3 |

