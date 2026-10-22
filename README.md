# SPM-OP

Parameter identification for a single particle model (SPM) of a lithium-ion cell, using particle swarm optimization (PSO) or a genetic algorithm (GA). Given current/voltage data from a "degraded" target cell, the code searches for the electrochemical parameters (and each profile's initial state) that minimize the voltage RMSE of the SPM.

This is a baseline for physics-informed neural network (PINN) parameter identification. It uses the SPM with the LiCoO2 parameter set from Scott Moura's SPMeT, not the LG M50 SPMe used in PINN-SPMe.

For model details, see Scott Moura's [SPMe observer paper](https://ecal.studentorg.berkeley.edu/pubs/SPMe-Obs-Journal-Final.pdf) and the [SPMeT repository](https://github.com/scott-moura/SPMeT). The parameter sets are based on SPMeT's [`params_LCO.m`](https://github.com/scott-moura/SPMeT/blob/master/param/params_LCO.m). The datasets were simulated with a modified SPMeT [`spme.m`](https://github.com/scott-moura/SPMeT/blob/master/spme.m) with the electrolyte dynamics disabled.

## Problem setup

- **Model:** SPM with finite-difference radial discretization (`Nr_model = 20`). The solid diffusion is linear and simulated with `lsim` on the state-space matrices from `spm_plant_obs_mats_paramID_{n,p}.m`. The electrolyte concentration is constant, and the voltage comes from `nonlinearSPMOutputVoltage_paramID.m`.
- **Known and target cells:**
  - `params_setting.m` is the known (fresh) cell.
  - `params_setting_2.m` is the target (degraded) cell: n_Li_s, k_n, k_p, D_s_n, D_s_p and D_e are scaled by 0.8, and R_f_n by 1.2.
  - Data for each cell are SPM simulations at `Nr = 100` (`data/SPM_FDM_nr100_simulation_{profile}_cell{1|2}.csv`).
- **Decision variables:** a subset `p_targets` of six parameters, each scaled by its true value: 1 n_Li_s, 2 R_f_n, 3 k_n, 4 k_p, 5 D_s_n, 6 D_s_p. Unselected parameters are fixed at the truth. Search bounds are 1/`bound` to `bound` times the known cell's value (`bound = 2`).
- **Initial state:** the initial cathode concentration of each profile is an extra unknown, unless `tell_csp0 = true`. The anode concentration follows from n_Li_s through `init_cs_paramID.m`.
- **Objective:** voltage RMSE against data with 3% Gaussian noise added to current and voltage.
- **Optimizers:** MATLAB `particleswarm` (PSO) or `ga` (GA), both with `UseParallel` enabled.

## Usage

Requires MATLAB with the Global Optimization and Parallel Computing toolboxes. Create the output folder once with `mkdir results`.

```matlab
% Joint fit over five profiles (UDDSx2, FUDS, BJDST, Charge, Charge2); method set at the top (1 GA, 2 PSO)
optimize_spm_multi_profiles      % starts parpool, writes results/x_est.csv

% Single-profile fits, looping over profiles and methods
paramID_wrapper                  % calls optimize_spm.m, writes result_<targets>.txt
```

Edit the configuration block at the top of either script to change these settings:

- `cell_target`, `cell_known`
- `bound`, `profiles`, `p_targets`, `tell_csp0`
- `Nr_model`, `Nr_data`, `method`

On a SLURM cluster, `matlab_batch.sh` runs `optimize_spm_multi_profiles`. Edit the module name and paths for your cluster first.

## Reference result

`paramID_results.xlsx` holds the joint five-profile identification of the degraded cell, with targets n_Li_s, k_n, k_p, D_s_n and D_s_p. Values are estimate / true value (1.0 = exact), and R_f_n was fixed. The PINN row is the estimate recorded in the same sheet for the companion PINN-SPM model.

| Method | n_Li_s | k_n | k_p | D_s_n | D_s_p |
|---|---|---|---|---|---|
| PSO | 1.156 | 1.431 | 0.629 | 0.625 | 1.914 |
| GA | 1.253 | 0.746 | 0.826 | 0.628 | 1.971 |
| PINN | 0.985 | 1.086 | 0.707 | 2.396 | 1.784 |

The sheet also lists the identified initial-state ratios of each profile for PSO and GA, along with the relative errors.

## Files

| File | Purpose |
|---|---|
| `optimize_spm_multi_profiles.m` | Script: joint identification over several profiles |
| `optimize_spm.m` | Function: identification on one profile (used by the wrapper) |
| `paramID_wrapper.m` | Loops over profiles and methods and tabulates the results |
| `params_setting.m`, `params_setting_2.m` | Known (fresh) and target (degraded) parameter sets |
| `spm_plant_obs_mats_paramID_n.m`, `spm_plant_obs_mats_paramID_p.m` | Solid-diffusion state-space models per electrode |
| `nonlinearSPMOutputVoltage_paramID.m` | Voltage for a time series |
| `init_cs_paramID.m` | Initial concentrations from the candidate n_Li_s and c_s,p0 |
| `exch_cur_dens.m`, `refPotentialAnode.m`, `refPotentialCathode.m` | Kinetics and open-circuit potentials |
| `matlab_batch.sh` | SLURM job script |
| `data/` | SPM simulations at Nr = 20 and Nr = 100 for cells 1 and 2 (no header; columns: time, current, voltage, surface/average anode and cathode concentration) |

## Related repositories

- [SPM-EKF](https://github.com/yuichikajiura/SPM-EKF): EKF state estimation for the same SPM. It uses an n_Li_s value identified here as a pseudo-measurement.
- [PINN-SPMe](https://github.com/yuichikajiura/PINN-SPMe) and [LGM50](https://github.com/ESCL-OU/LGM50): PINN state and SEI estimation for an SPMe with degradation, and the MATLAB simulator that generates its data.

## Credits

The model functions are adapted from Scott Moura's SPMeT, which is licensed under GPL-3.0.
