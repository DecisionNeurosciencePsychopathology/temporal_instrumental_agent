# ![sceptic]()

# StrategiC Exploration/Exploitation of a Temporal Instrumental Contingency

This repository contains code realted to `SCEPTIC`, which was developed to model how humans explore and learn a complex reinforcement-based timing task. For more information see Dombrovski Hallquist [PAPER Selective Maintenance and Entropy-Driven Exploration Facilitate Human Learning of Temporal Instrumental Contingencies]().

## Requirements
Please note Matlab is required to run SCEPTIC 

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

## Team

[![Alex Dombrovski]()]() | [![Michael Hallquist]()]() | [![Jonathan Wilson]()](https://github.com/wilsonj3)
---|---|---
[Alex Dombrovski](https://github.com/dombrovski) | [Michael Hallquist](https://github.com/michaelhallquist) | [Jonathan Wilson](https://github.com/wilsonj3)


## License

