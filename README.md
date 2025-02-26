# SPM-OP

## Single Particle Model(SPM) Optimization and Parameterization

### Method
This is a set of code to identify parameters of SPM from battery cycling datasets using metaheuristic algorithm (Genetic Algorithm or Particle Swarm Optimization).
SPM is discretized by finite difference method is used as the battery model. For details, refer to [this paper](https://ecal.studentorg.berkeley.edu/pubs/SPMe-Obs-Journal-Final.pdf) from Professor Scott Moura at UC Berkeley, as well as [this repository](https://github.com/scott-moura/SPMeT). (Note: the electrolyte dynamics is neglected/disabled for the moment)

### Features
- Six model parameters (e.g. diffusion coefficient at the anode) can be set to be identified 
- Two cell models (fresh and aged) included, of which parameters based on [this](https://github.com/scott-moura/SPMeT/blob/master/param/params_LCO.m)
- Several cycling datasets are prepared by simulating SPM as coded [here](https://github.com/scott-moura/SPMeT/blob/master/spme.m) (with electrolyte dynamics disabled and some modifications)

### Usage
To identify parameter from any one of current profile, run:
```console
matlab param_ID_wrapper.m
```
To choose the parameter to be identified or limit the range of parameter search (as a ratio to the parameters of known cell), or more, you can change the configuration in the first several lines of the codes.

To identify parameter using multiple current profiles at once, run:
```console
matlab optimize_spm_multi_profiles.m
```
Again, you can change configuration of the code to choose parameter to be identified or more. 
