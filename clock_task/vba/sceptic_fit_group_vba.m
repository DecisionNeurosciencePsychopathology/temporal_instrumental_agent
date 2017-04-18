function sceptic_fit_group_vba(sceptic_config_file)
%Parent script for VBA inversion of SCEPTIC models for group level analysis and processing
%Requirements:
%Must have VBA toolbox (J. Daunizeau, V. Adam, L. Rigoux) downloaded and installed 
%see here for more help: https://mbb-team.github.io/VBA-toolbox/
%
%Must create configuration file with "create_scecptic_configuration_struct" script

clearvars -except sceptic_config_file
%curpath = fileparts(mfilename('fullpath'));

%Load in configuration file
try
    s=load(sceptic_config_file);
    tmp=fieldnames(s);
    s = s.(tmp{:});
catch
    error('Configuration file was not able to be loaded')
end

%Is VBA toolbox on path and/or downloaded
try which VBA_NLStateSpaceModel; catch, error('VBA toolbox is not on path and/or... downloaded, see script description for more info.'); end

%% Set up file paths
behavfiles = glob(s.behavfiles); %Path of subjects data
group_dir = s.group_dir; %Path to save group posterior results 

%% chose models to fit
%Models to choose from are:
%1) fixed
%2) fixed_decay
%3) fixed_uv
%4) kalman_softmax
%5) kalman_uv_sum
modelnames = s.modelnames;

% get ID list
id_list = NaN(length(behavfiles),1);

%% main loop
%Initialize posterior log evidence array
L = NaN(length(modelnames),length(behavfiles));

%Start main model loop
for m=1:length(modelnames)
    model = char(modelnames(m));
    parfor sub=1:length(behavfiles)
        
        %Pull subject's data file path
        data_file = behavfiles{sub};
        
        %Pull subject id from file name
        id=regexp(data_file,'\d{5,6}','match');
        id=str2double(id{:});
        fprintf('Fitting %s subject %d \r',model,id)
        
        %Run model through VBA toolbox
        [posterior,out] = clock_sceptic_vba(s,id,model,data_file);
        L(m,sub) = out.F;
        id_list(sub) = id;
    end
end
cd(group_dir);
filename = sprintf('SHIFTED_U_grp_L_%d_nbasis%d_nsteps%d_uaversion_not_allModels_fixed_prop_spread_frank_all_subjs',s.nbasis,s.n_steps,s.u_aversion);
save(filename, 'L', 'id_list');

