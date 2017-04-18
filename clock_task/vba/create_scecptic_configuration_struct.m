function s=create_scecptic_configuration_struct(save_flag)
%Configuration file for running parent script (sceptic_group_fit_vba.m), which is a group level analysis script using the SCEPTIC model and the VBA toolbox

%If you want to save the conifguration file (saved by default)
if nargin<1; save_flag=1; end

%TASK NAME
s.task_name = 'hallquist_clock'; %More for saving results purposes

%FILE I/O
%behavfiles  = Path to subject data
%results_dir = Path to individual VBA results (if you want to save them)
%group_dir   = Path to group posteriors (if you want to save them)

s.behavfiles='';
s.results_dir='';
s.group_dir='';

%MODELS
%Models to choose from are:
%1) fixed
%2) fixed_decay
%3) fixed_uv
%4) kalman_softmax
%5) kalman_uv_sum
%model names must be wrapped in qoutes (i.e. 'fixed' 'fixed_decay' ect)

s.modelnames = {'fixed_decay'};

%MODEL PARAMETERS
%nbasis                   = Number of basis functions                                                                               | Default = 24
%multinomial              = VBA option defning the likelihood function of what you are tyring to predict                            | Default = 1
%multisession             = VBA option if you would like to vary/keep constant parameters during muliple ecxperimental sessions     | Default = 0
%fixed_params_across_runs = If when using multisession you would like to keep all model parameters fixed across multiple runs       | Default = 0
%fit_propspread           = If you would like the toolbox to fit the prop_spread parameter                                          | Default = 0
%n_steps                  = Number of time bins                                                                                     | Default = 50
%u_aversion               = Allow for uncertainty aversion in UV_sum                                                                | Default = 1
%saveresults              = If you would like to save the results                                                                   | Default = 0
%graphics                 = If you would like to display graphics                                                                   | Default = 0
%cstruct                  = Take std over all timesteps and possible draws per condition for sigma_noise parameter                  | Default = []
%range_RT                 = Maximum time point of reaction time range (on bin scale (ex 4000ms in 10's time bin scale is 400))      | Default = 400

s.nbasis=16;
s.multinomial=1;
s.multisession=0;
s.fixed_params_across_runs=0;
s.fit_propspread=0;
s.n_steps=50;
s.u_aversion=1;
s.saveresults=1;
s.range_RT=500;

%SAVING
s.results_str = 'clock_task_%s_%d_%s'; %Args are is task_name, id, model
s.results_str_values =  s.task_name;
if save_flag
   [filename, pathname] = uiputfile('sceptic_data_struct.mat','Save file name');
   save([pathname filename], 's');
end