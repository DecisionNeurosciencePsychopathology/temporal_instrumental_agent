function [p_sub,o_sub,p_group,o_group] = run_clock_vba_mfx(vba_df)



%Parse out the vba dataframe
y = vba_df.y;
u = vba_df.u;
%options = vba_df.options;
%dim = options{1}.dim;

%models = {'kalman_logistic', 'kalman_uv_sum', 'kalman_softmax'};
% models = {'fixed', 'fixed_decay', 'fixed_uv', 'kalman_softmax', 'kalman_logistic', 'kalman_uv_sum'};
models = {'fixed_decay'};


%VBA mfx does not like unequal trials, bring this up to Alex
%Find which index have less than max trials
trial_lengths = cellfun(@length,vba_df.y);
missing_data_index = find(max(trial_lengths) > trial_lengths);

%Remove missing subjects 
y(missing_data_index) = [];
u(missing_data_index) = [];

for model = models
    
    close all; %Get rid of extensive gpu memory figs
    
    %Initialize the output matrix
    vba_mfx_df = table();
    
    %Load in the options 
    %options = load([unix_str('E:\data\sceptic\vba_out\specc_clock\vba_mfx_input\') sprintf('vba_mfx_input_%s_specc',model{:})]);
%     options = load([unix_str('E:\data\sceptic\vba_out\explore_clock\vba_mfx_input\') sprintf('vba_mfx_input_%s_explore',model{:})]);
    options = vba_df.options;
    dim = options{1}.dim;
    
    options(missing_data_index) = [];
    
    %Designate the f and g function handles
    f_name = options{1}.f_fname; %Evolution function
    g_name = options{1}.g_fname; %Observation function

    %clear vars for new output
    clearvars p_sub o_sub p_group o_group
    
    %[p_sub,o_sub,p_group,o_group] = VBA_MFX(y,u,f_fname,g_fname,dim,options,priors_group);
    [p_sub,o_sub,p_group,o_group] = VBA_MFX(y,u,f_name,g_name,dim,options);
    
    %Grab the outputs
%     vba_mfx_df.ID = vba_df.ID;
%     vba_mfx_df.p_sub = p_sub;
%     vba_mfx_df.o_sub = o_sub;
%     
%     %Try to save as one data frame but break it up just in case - make sure
%     %path is right!
%     save(sprintf('E:/data/sceptic/vba_out/explore_clock/vba_mfx_input/vba_mfx_output_vba_mfx_df_%s',model{:}),'vba_mfx_df', '-v7.3')
    save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/n76/vba_mfx_input/vba_mfx_output_p_sub_%s',model{:}),'p_sub', '-v7.3')
    save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/n76/vba_mfx_input/vba_mfx_output_o_sub_%s',model{:}),'o_sub', '-v7.3')
    save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/n76/vba_mfx_input/vba_mfx_output_p_group_%s',model{:}),'p_group', '-v7.3')
    save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/n76/vba_mfx_input/vba_mfx_output_o_group_%s',model{:}),'o_group', '-v7.3')
    
%     save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/explore_clock/vba_mfx_input/vba_mfx_output_p_sub_%s',model{:}),'p_sub', '-v7.3')
%     save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/explore_clock/vba_mfx_input/vba_mfx_output_o_sub_%s',model{:}),'o_sub', '-v7.3')
%     save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/explore_clock/vba_mfx_input/vba_mfx_output_p_group_%s',model{:}),'p_group', '-v7.3')
%     save(sprintf('/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/explore_clock/vba_mfx_input/vba_mfx_output_o_group_%s',model{:}),'o_group', '-v7.3')
%     
end