%loads in subjects' data and fits SCEPTIC models using VBA;
close all; clear;
    
% Note: this function looks for 'sceptic_dataset' and 'sceptic_model'
% as environment variables so that this script can be scaled easily for batch processing

so = []; %you can set custom sceptic options here that will override the function below.

%for manual tests
%so.model = 'fixed_uv';
%so.dataset = 'mmclock_fmri';
%so.model='decay';
%so.dataset='explore';
%so.sceptic_dataset=so.dataset;
%setenv('matlab_cpus', '2');

so = sceptic_validate_options(so); %initialize and validate sceptic fitting settings
so.mfx = 1; %for mfx subdir naming

%setup paths, parallel pools, etc. based on user environment
[so, poolobj, behavfiles] = sceptic_setup_environment(so);

%extract IDs for record keeping

if(strcmpi(so.dataset,'explore'))
  ids=cellfun(@getexploreid,behavfiles,'UniformOutput', false);
else
  [~,fnames]=cellfun(@fileparts, behavfiles, 'UniformOutput', false);
  ids=cellfun(@(x) char(regexp(x,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match')), fnames, 'UniformOutput', false);
end

%zero pad IDs if lengths vary so that character strings have consistent length for concatenation of group statistics
maxlen = max(cellfun(@(x) strlength(x), ids));

ids=cellfun(@(x) sprintf('%0*s', maxlen, x), ids, 'UniformOutput', false);

%puts a nested cell in each element
%ids=regexp(fnames,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

ns = length(behavfiles);
y_all = cell(ns, 1);
u_all = cell(ns, 1);
options_all = cell(ns, 1);

n_t=NaN(1,ns);

global rew_rng_state 
rew_rng_seed = 99;

for sub = 1:ns
    fprintf('Loading subject %d id: %s \n', sub, ids{sub});
    [data, y, u] = sceptic_get_data(behavfiles{sub}, so);
    [options, dim] = sceptic_get_vba_options(data, so);
    n_t(sub) = dim.n_t; % allow for variation in number of trials across subjects
    
    %Set up sigma noise for every point in u or hidden state?
    rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
    rew_rng_state=rng;
    [~,idx] = unique(data.run);
    conditions=data.rewFunc(idx);
    sigma_noise = zeros(length(conditions), 1);
    run_length = dim.n_t/options.inF.n_runs;
    for i = 1:length(conditions)
      sni = repmat(std(arrayfun(@(x) RewFunction(x*100, conditions(i), 0), options.inF.tvec))^2, 1, run_length);
      sigma_noise(i) = mean(sni);
    end
    sigma_noise = mean(sigma_noise);
    options.inF.sigma_noise = sigma_noise;
    
    %populate sigmas with sigma_noise as prior
    if contains(so.model, 'fixed_uv')
      hidden_state_index=1:so.hidden_states*so.nbasis; %total number of hidden states (inF.hidden_states is the number of state vectors)
      hidden_state_index = reshape(hidden_state_index, so.nbasis, so.hidden_states); %3 x nbasis here
      
      %uncertainty is second hidden state
      options.priors.muX0(hidden_state_index(:,2)) = sigma_noise;
    end
    
    % populate data structures for VBA_MFX
    y_all{sub} = y;
    u_all{sub} = u;
    options_all{sub} = options;
end

% add the n_t vector to dim before passing to MFX
dim.n_t=n_t;

%options for MFX
options_group.TolFun=1e-2;
options_group.MaxIter=50;
options_group.DisplayWin=0;
options_group.verbose=1;

priors_group.muPhi = zeros(dim.n_phi,1); %temperature -- exp(phi(1))
priors_group.SigmaPhi = 1e1*eye(dim.n_phi); %variance on temperature (before exponential transform)
priors_group.muTheta = zeros(dim.n_theta,1); %learning rate (alpha), selective maintenance (gamma) -- before logistic transform
priors_group.SigmaTheta =  1e1*eye(dim.n_theta); %variance of 10 on alpha and gamma

priors_group.muX0 = zeros(so.nbasis*so.hidden_states,1); %have PE and decay as tag-along states
priors_group.SigmaX0 = zeros(so.nbasis*so.hidden_states, so.nbasis*so.hidden_states); %have PE and decay as tag-along states

%for initial hidden state values, need to populate sigmas for fixed_uv variants at *group* level
%in fact, MFX will overwrite the priors.muX0 with group values, undoing subject-wise setup above
if contains(so.model, 'fixed_uv')
  hidden_state_index=1:so.hidden_states*so.nbasis; %total number of hidden states (inF.hidden_states is the number of state vectors)
  hidden_state_index = reshape(hidden_state_index, so.nbasis, so.hidden_states); %3 x nbasis here
  
  %uncertainty is second hidden state
  priors_group.muX0(hidden_state_index(:,2)) = sigma_noise;
end

priors_group.a_vX0 = Inf([1, so.nbasis*so.hidden_states]); %use infinite precision prior on gamma for X0 to treat as fixed (a = Inf; b = 0)
priors_group.b_vX0 = zeros([1, so.nbasis*so.hidden_states]);

[p_sub, o_sub, p_group, o_group] = VBA_MFX_parallel(y_all, u_all, so.evo_fname, so.obs_fname, dim, options_all, priors_group, options_group);
%[p_sub, o_sub, p_group, o_group] = VBA_MFX(y_all, u_all, @h_sceptic_fixed_decay_fmri, @g_sceptic, dim, options_all, priors_group, options_group);

delete(poolobj);

%populate subject ids into output.options.inF structure since these are added to subject statistics below
for s=1:length(o_sub)
    o_sub{s}.options.inF.id = ids{s};
    o_sub{s}.options.inG.id = ids{s};
    
    %populate ffx parameters from o_group structure to p_sub structure for extraction
    p_sub{s}.ffx = o_group.initVBA.p_sub{s};
    p_sub{s}.ffx = rmfield(p_sub{s}.ffx, {'SigmaX', 'iQx', 'iQy'}); %large matrices not needed for anything (they use lots of disk space)
end

%create a structure with just the barebones useful parameters from subjects
s_all = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub, o_sub, 'UniformOutput', false);

%output MFX results

%save group outputs
save(sprintf('%s/group_fits_%s_%s', so.output_dir, so.model, so.dataset), 'ids', 'so', 's_all');

[group_global, group_trial_level] = extract_group_statistics(s_all, ...
    sprintf('%s/%s_%s_mfx_sceptic_global_statistics.csv', so.output_dir, so.dataset, so.model), ...
    sprintf('%s/%s_%s_mfx_sceptic_trial_outputs_by_timestep.csv', so.output_dir, so.dataset, so.model));

%save the basis as a csv file
vnames = cellfun(@(x) strcat('Time_', num2str(x)), num2cell(1:so.ntimesteps), 'UniformOutput', false);
basis_mat = array2table(o_sub{1}.options.inF.gaussmat, 'VariableNames', vnames);

writetable(basis_mat, sprintf('%s/%s_%s_mfx_sceptic_basis.csv', so.output_dir, so.dataset, so.model));

%save group outputs
save(sprintf('%s/group_fits_%s_%s', so.output_dir, so.model, so.dataset), 'ids', 'so', 's_all', 'group_global', 'group_trial_level');

%too huge to save into one .mat file
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_psub.mat'], 'p_sub', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_pgroup.mat'], 'p_group', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_ogroup.mat'], 'o_group', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_osub.mat'], 'o_sub', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_settings.mat'], 'priors_group', 'options_group', 'y_all', 'u_all', 'ids', '-v7.3');

%free energy (log evidence)
L = o_group.within_fit.F;

% save just the log evidence, L
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_L.mat'], 'L', '-v7.3');

