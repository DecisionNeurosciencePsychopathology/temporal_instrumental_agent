# ![sceptic]()

# StrategiC Exploration/Exploitation of a Temporal Instrumental Contingency

This repository contains code realted to `SCEPTIC`, which was developed to model how humans explore and learn a complex reinforcement-based timing task. For more information see Dombrovski Hallquist [Selective Maintenance and Entropy-Driven Exploration Facilitate Human Learning of Temporal Instrumental Contingencies]().

## Requirements
Please note Matlab is required to run `SCEPTIC`

You must have the VBA tool box downloaded and installed, for instructions see [here](https://mbb-team.github.io/VBA-toolbox/download/) for more information see [here](https://mbb-team.github.io/VBA-toolbox/wiki/)

## Install

Clone or download this repo

## Usage

1.) Create a config file by running:

```
s=create_scecptic_configuration_struct
```

2a.) Execute the parent function to run model(s) on all subjects
```
sceptic_fit_group_vba(sceptic_config_file)
```

2b.) Conversely run model(s) on a single subject
```
[posterior,out] = clock_sceptic_vba(s,id,model,data_file)
``` 

3.) Compare models

```matlab
sceptic_grp_BMC
```

## Current models

##### fixed

Fixed learning rate, value-based choice

##### fixed decay

selective maintenance, fixed learning rate, value-based choice, action values decay

##### fixed uv

fixed learning rate, choice controlled by weighted sum of value and uncertainty

##### kalman softmax

 Kalam filter (KF) learning rule, value-based choice

##### kalman uv sum

KF learning rule, choice controlled by weighted sum of value and uncertainty

## Configuration file

##### `task name`

Which version of the clock task is being analyzed Default: hallquist_clock

##### `behavfiles`

Path pointing to processed subjects' data Format: csv 

##### `results_dir`

Path pointing to where to save indivdual results outputted from VBA toolbox Outputs: `posterior` `out`

##### `group_dir`

Path pointing to where to save the group results 

##### `modelnames`

Cell containing strings of model names Default: `{'fixed' 'fixed_decay' 'fixed_uv' 'kalman_softmax' 'kalman_uv_sum'}`

##### `nbasis`

Number of basis functions Default: `24`

##### `multinomial`

VBA option defning the likelihood function of what you are tyring to predict Default: `1`

##### `multisession`

VBA option if you would like to vary/keep constant parameters during muliple ecxperimental sessions Default: `0`

##### `fixed_params_across_runs`

If when using multisession you would like to keep all model parameters fixed across multiple runs Default `0`

##### `fit_propspread`

If you would like the toolbox to fit the prop_spread parameter Default: `0`

##### `n_steps`

Number of time bins Default:`50`

##### `u_aversion`

Allow for uncertainty aversion in UV_sum Default: `1`

##### `saveresults`

If you would like to save the results Default: `0`

##### `graphics`

If you would like to display graphics `0`

##### `cstruct`

Take std over all timesteps and possible draws per condition for sigma_noise parameter Default: `[]`

##### `range_RT`

Maximum time point of reaction time range (on bin scale (ex 4000ms in 10's time bin scale is 400)) Default: `400`

## Team

[![Alex Dombrovski]()]() | [![Michael Hallquist]()]() | [![Jonathan Wilson]()](https://github.com/wilsonj3)
---|---|---
[Alex Dombrovski](https://github.com/dombrovski) | [Michael Hallquist](https://github.com/michaelhallquist) | [Jonathan Wilson](https://github.com/wilsonj3)


## License

