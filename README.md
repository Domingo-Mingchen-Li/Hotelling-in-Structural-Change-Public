# Hotelling in Structural Change: Green Paradox and Delayed Transformation

[![MATLAB](https://img.shields.io/badge/MATLAB-R2023a%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Author:** Mingchen Li  

## Overview
This repository contains the replication codes for the working paper *"Hotelling in Structural Change: Green Paradox and Delayed Transformation"*. 

## Repository Structure

### 1. Core Solvers
* `solve_shooting_2sector_inverted.m`: The core algorithm featuring inverted physical logic and smart damping for solving non-linear transitional dynamics.
* `run_precision_homotopy_2sector.m`: An adaptive multi-dimensional homotopy engine for policy interventions and parameter transitions.
* `solve_backward_step_2sector.m`: The single-period general equilibrium solver.
* `solve_full_backward_path_2sector.m`: Executes the backward-shooting process to derive the entire transitional dynamics path given the terminal T-period equilibrium configuration.
* `solve_terminal_ABGP_2sector.m`: Computes the terminal Asymptotic Balanced Growth Path (ABGP).
* `static_kernel_2sector.m`: The static equilibrium kernel resolving factor prices and allocations on ABGP.
* `get_growth_rates_2sector.m`: Calculates the ABGP growth rates.

### 2. Experiments & Tests
* `test_policy_intervention.m`: Simulates the implementation of intuitively clean industrial policies and visualizes the resulting "Structural Green Paradox".
* `test_counterfactual.m`: A counterfactual experiment completely freezing the structural change mechanism to isolate the pure Baumol/relative price effect from the Engel/income effect.

### 3. Data & Visualization
* `plot_Figure1_Baseline.m`: Renders the baseline macroeconomic transitional dynamics.
* `plot_Figure2_PolicyDeviations.m`: Visualizes the deviation of macroeconomic variables under policy shocks.
* `analyze_green_paradox_decomposition.m`: Performs a novel 3-way log-difference decomposition of the Green Paradox.

---

## How to Run & Replicate Results

### 1. Running Custom Simulations
You can generate your own custom policy scenarios by modifying and executing the two test scripts (`test_policy_intervention.m` and `test_counterfactual.m`) located in the main program directory. 

*Note: Due to the highly non-linear nature of the dynamic general equilibrium model, not all arbitrary policy combinations are robust or guaranteed to converge. Successfully simulated results will be automatically saved to the `Sim_Results/` directory, which can subsequently be visualized using the provided plotting scripts.*

### 2. Replicating Paper Results
To facilitate exact replication of the results and figures presented in the paper, pre-computed simulation data has been provided in the `Sim_Results/` folder. You can replicate the specific figures and quantitative analyses by executing the following scripts:

* **Figure 1 (Baseline Dynamics):** Run `plot_Figure1_Baseline.m`. You can input the path to *any* generated `.mat` data file within the `Sim_Results/` folder to render the baseline transition paths.
* **Figure 2 (Policy Deviations):** Run `plot_Figure2_PolicyDeviations.m` and ensure the file path within the script is set to `Sim_Results/service_subsidy_10_periods.mat`.
* **Green Paradox Decomposition:** Run `analyze_green_paradox_decomposition.m`. Inputting `Sim_Results/service_subsidy_10_periods.mat` (or other policy shock data files) will accurately replicate the exact log-difference decomposition results highlighted in the paper.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
