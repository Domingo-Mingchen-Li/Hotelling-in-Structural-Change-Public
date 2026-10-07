# The Green Paradox in Structural Change

MATLAB code and saved results for the paper's transition, calibration, robustness, counterfactual and welfare exercises. Requires MATLAB R2024a or later.

## Reproduce the paper

Set the repository root as the MATLAB current directory and run:

```matlab
reproduce_main_figures
reproduce_appendix_figures
```

The first command writes Figure1–Figure16 and Table2 to `outputs/paper_figures/`. The second writes FigureA1–FigureA4 and TableA2–TableA7 to `outputs/appendix_figures/`. Figures are exported as PNG and vector PDF; tables as CSV and LaTeX. Both commands use the bundled results in `results/paper/`.

## Figure and code correspondence

| Paper figures | Section | Exporter | Contents |
|---|---|---|---|
| Figure1 | 5.1 | `export_paper_figures.m` | Baseline transition |
| Figure2–4 | 5.2 | `export_paper_figures.m` | Temporary services subsidy |
| Figure5 | 5.3 | `export_paper_figures.m` | Permanent services subsidy |
| Figure6–8 | 5.4 | `export_paper_figures.m` | Consumption tax and investment subsidies |
| Figure9 | 5.5 | `export_paper_figures.m` | Resource tax and labour allocation |
| Figure10 | 5.6 | `export_paper_figures.m` | Cumulative extraction over dates 0–40 |
| Figure11–13 | 5.7 | `export_counterfactual_writing.m` | Flat resource-input shares: two policy responses and baselines |
| Figure14–16 | 5.8 | `export_policy_evaluation_writing.m` | Industrial-policy welfare, optimal resource tax and shared-history comparison |
| FigureA1–2 | Appendix D | `export_paper_figures.m` | Horizon comparisons |
| FigureA3–4 | Appendix E | `export_measurement_writing.m` | Measurement robustness |

| Paper table | Contents | Exporter |
|---|---|---|
| Table2 | Cumulative extraction over 10, 50, 100 and 200 flows | `export_paper_cumulative_table.m` |
| TableA2–A3 | Recalibrated parameters and ten-flow cumulative extraction | `export_measurement_writing.m` |
| TableA4–A7 | Persistence kernel, damage scenarios, industrial-policy welfare and optimal-tax comparison | `export_policy_evaluation_writing.m` |

Edit `config/paper_figure_config.m` to select result files, output folders and display windows. Layout and labels are in the exporters. Table1 and TableA1 report the fixed calibration and industry mapping given in the paper and input files.

## Solver and paper correspondence

`core/transition_static_kernel.m` implements the within-period equilibrium in Section 3. `transition_intertemporal_residual.m` implements the Euler and Hotelling conditions. The global residual, terminal boundary and sparse Newton routines correspond to Section 4.1 and Appendix B.

The previous version of the paper used backward shooting. Repeated multiplication of period-by-period Jacobians amplified sensitivities along the transition path and made the method numerically unstable. This revision replaces backward shooting with a sparse solver that solves the entire path simultaneously, without relying on these Jacobian products. The previous solver is available on the [backward-shooting-v1 branch](https://github.com/Domingo-Mingchen-Li/Hotelling-in-Structural-Change-Public/tree/backward-shooting-v1). The current solver jointly determines `z_t=(K_t,R_t,r_t,h_t,E_t)` over dates 0,…,T, with other prices and allocations reconstructed within each period. Five equations per preterminal date, two initial-stock conditions and three ABGP terminal conditions close the system. Logarithmic coordinates, a sparse Jacobian and damped Newton steps handle the full path. Separate numerical audits assess derivatives, equilibrium conditions, restarts and horizon sensitivity.

The paper uses unscaled technology growth factors and a revised exact three-moment calibration. The retained `legacy` equation convention uses the same sparse numerical architecture.

## Recompute the results

Configurations are in `config/`; fresh runs are saved in `outputs/`. The initial guess is `data/paper_initial_paths.mat`.

```matlab
setup_solver
run_validation
run_experiments
run_robustness
report = run_calibration();
```

`solver_config.m` and `policy_config.m` specify the transition parameters and policies. After recalibration, transfer the fitted parameters and warm-start path to the experimental configuration.

Measurement robustness (Appendix E):

```matlab
cal = run_measurement_robustness_calibration();
report = run_measurement_robustness();
cal = run_measurement_inputs_calibration();
report = run_measurement_inputs();
cal = run_measurement_mapping_expanded_calibration();
report = run_measurement_mapping_expanded();
```

These implement annual-average shares, domestic inputs and the C26 mapping. The mapping workflow includes its continuation seed and calibration bounds.

Resource-input heterogeneity (Section 5.7):

```matlab
report = run_counterfactual();
```

Climate and welfare (Section 5.8 and Appendix F), in execution order:

```matlab
climate = run_climate_parameterization();
welfare = run_industrial_welfare();
ramsey = run_ramsey();
shared = run_ramsey_shared_state_counterfactual();
```

The shared-history experiment switches input coefficients at L=10. Each stage uses the preceding stage's saved handoff. Climate and optimal-policy workflows include T=300 checks. Default runs use compact logging; numerical options remain editable in the configs.

## Conventions and data

Temporary transition policies apply at dates 0–9. Welfare exercises use an unexpected announcement at L=10 and perfect foresight thereafter. A cumulative window of H flows covers dates 0,…,H−1. Policy responses use each economy's own baseline. The timing decomposition reports an equilibrium accounting identity in log points.

The calculations use a finite-horizon ABGP approximation. Environmental welfare is a normative scenario in normalized model units; Ramsey results are local solutions under the configured control bounds and terminal policy. The Pigouvian curve evaluates the analytical condition on the numerically optimized allocation.

`data/measurement/` holds China aggregates and industry mappings. `data_generation/` contains the Python scripts to rebuild them. With the original WIOT and SEA workbooks available:

```bash
python data_generation/reproduce_beta.py --input-dir /path/to/WIOT --output-dir /path/to/derived
python data_generation/build_measurement_parameters.py --sea /path/to/SEA.xlsx --wiot-dir /path/to/WIOT --output-dir /path/to/derived
python data_generation/build_sector_mapping_parameters.py --data-dir data/measurement --output-dir /path/to/derived_mapping
```

The SEA reader requires `openpyxl`. WIOT filenames run from `WIOT2000_Nov16_ROW.xlsb` through `WIOT2014_Nov16_ROW.xlsb`; SEA uses a `DATA` sheet with country, variable, industry code and year columns. `results/PROVENANCE.json` records the saved result sources and checksums.

Result integrity, dependencies and MATLAB syntax have been checked. Figure rendering remains to be checked in MATLAB. `.gitignore` excludes generated outputs.
