# ![sceptic]()

#SCEPTIC

Description

## Requirements
These scripts run in Matlab

You must have the VBA tool box downloaded and installed, for instructions see [here](https://mbb-team.github.io/VBA-toolbox/download/) for more information see [here](https://mbb-team.github.io/VBA-toolbox/wiki/)

## Install

Clone or download this repo

## Usage

1.) Create a config file by running:

```matlab
s=create_scecptic_configuration_struct
```

2a.) Execute the parent function to run model(s) on all data
```matlab
sceptic_fit_group_vba(sceptic_config_file)
```

2b.) Conversely run model(s) on a single data point
```matlab
[posterior,out] = clock_sceptic_vba(s,id,model,data_file)
``` 

3.) Compare models

```matlab
sceptic_grp_BMC
```

## Current models

##### fixed

Description

##### fixed decay

Description

##### fixed uv

Description

##### kalman softmax

Description

##### kalman uv sum

Description

## Team

[![Alex Dombrovski]()]() | [![Michael Hallquist]()]() | [![Jonathan Wilson]()](https://github.com/wilsonj3)
---|---|---
[Alex Dombrovski](https://github.com/dombrovski) | [Michael Hallquist](https://github.com/michaelhallquist) | [Jonathan Wilson](https://github.com/wilsonj3)


## License

