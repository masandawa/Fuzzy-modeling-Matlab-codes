# Intensity-dependent fuzzy modelling of malaria–hookworm co-infection

MATLAB code accompanying the manuscript

> Masandawa L, Mbalawata IS, Kinung'hi S, Mirau SS. *Intensity-dependent fuzzy modelling of malaria–hookworm co-infection to assess integrated deworming and malaria control towards the WHO 2030 targets.*

The code reproduces all numerical results and the seven data figures of the manuscript.

## Quick start

1. Download or clone this repository and open its folder in MATLAB.
2. Open `run_all.m` and press **Run**.

The key results are printed in the Command Window. The seven figures then open on screen and are saved as PDF files in `figures/`.

Settings at the top of `run_all.m`:

| Setting | Effect |
|---|---|
| `RUN_ANALYSIS = true` | Runs the full analysis (about 30–60 minutes), then draws the figures |
| `RUN_ANALYSIS = false` | Draws the figures in seconds from the saved results in `results/figure_data.mat` |
| `WRITE_CSV = true` | Also saves CSV tables of all results in `results/` |
| `RUN_CHECKS = true` | Also runs the robustness checks (several hours) |

**Requirements:** MATLAB with the Optimization Toolbox (for `fsolve`). The code also runs in GNU Octave 8 or later.

## Files

| File or folder | Content |
|---|---|
| `run_all.m` | **Run this file.** Runs the analysis and draws the figures |
| `main_fuzzy_malaria_hookworm.m` | The analysis, covering:<br>• calibration and model evaluation<br>• threshold quantities, fuzzy reproduction numbers and identifiability<br>• bifurcation and global convergence<br>• critical deworming coverage<br>• sensitivity analysis<br>• fuzzy evaluation of six control scenarios for 2026–2030 |
| `make_figures.m` | Draws the seven data figures of the manuscript |
| `src/` | Model functions used by the analysis |
| `checks/` | Optional robustness checks (Section 5.5 of the manuscript) |
| `results/figure_data.mat` | Results used for the figures in the manuscript |

### Figures produced

| Manuscript figure | Content | Saved as |
|---|---|---|
| Fig. 1 | Fuzzy intensity groups and intensity-dependent parameters | `fig2_membership.pdf` |
| Fig. 2 | Model flow diagram (drawn in the LaTeX manuscript) | – |
| Fig. 3 | Model evaluation | `fig6_validation.pdf` |
| Fig. 4 | Bifurcation diagram and global convergence | `fig_bifurcation.pdf` |
| Fig. 5 | Reproduction numbers versus intensity | `fig3_R0.pdf` |
| Fig. 6 | Fuzzy equilibria | `fig4_fuzzy_equilibria.pdf` |
| Fig. 7 | Critical deworming coverage | `fig_coverage.pdf` |
| Fig. 8 | Fuzzy scenario trajectories, 2026–2030 | `fig8_scenarios.pdf` |

### Model functions (`src/`)

| Group | Functions |
|---|---|
| Model and parameters | `fmh_params` (baseline parameters), `fmh_shape`, `fmh_resp` (intensity-response functions), `fmh_derived`, `fmh_rhs` (model equations) |
| Thresholds and equilibria | `fmh_R0`, `fmh_invasion` (invasion factor and invasion number), `fmh_boundary`, `fmh_equilibrium`, `fmh_eqrobust`, `fmh_ic`, `fmh_jac`, `fmh_red`, `fmh_maxeig`, `fmh_bifcoef` (bifurcation coefficients) |
| Deworming thresholds | `fmh_critcov` (averaged treatment), `fmh_critcov_imp` (exact, periodic rounds) |
| Calibration and evaluation | `fmh_calibrate`, `fmh_hindcast`, `fmh_hindcast_at` |
| Simulation | `fmh_simulate` (periodic deworming rounds), `fmh_schedule`, `fmh_summaries`, `fmh_mhi` |
| Fuzzy tools | `fmh_trimf`, `fmh_alphacut`, `fmh_fev` (fuzzy expected value), `fmh_cutvals`, `fmh_fuzzyinputs`, `fmh_sample`, `fmh_design`, `fmh_static_out` |
| Sensitivity and output | `fmh_prcc`, `fmh_writecsv` |

## Reproducibility note

All random sampling uses fixed seeds. The saved results were computed with GNU Octave 8.4.

MATLAB and Octave use different random-number generators, so some numbers differ slightly between them:
- **Identical on both:** results that involve no random sampling, such as the calibration, reproduction numbers, critical coverage and core scenario outcomes.
- **Differ slightly:** sampled quantities, such as the bounds at α = 0, possibility/necessity values and sensitivity coefficients. The differences are within the variation documented by the robustness checks.

## Data

All data are published and listed in the manuscript (Table 4).
