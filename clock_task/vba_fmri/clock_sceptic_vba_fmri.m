function [posterior,out] = clock_sceptic_vba_fmri(data_file,model,n_basis,multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,u_aversion, saveresults, graphics)
%% fits SCEPTIC model to single-subject Clock Task subject data using VBA toolbox
% example call:
% [posterior,out]=clock_sceptic_vba(10638,'modelname',nbasis,multinomial,multisession,fixed_params_across_runs,fit_propsrpead)
% data_file:    CSV file containing raw trial-level data from clock task
% only works with 'fixed' (fixed learning rate SCEPTIC) so far
% n_basis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
% fixed_params_across_runs -- self-explanatory
% fit_propspread -- makes temporal generalization within the eligibility trace a free parameter
% n_steps:      number of time bins
% u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum

close all

%% uncertainty aversion for UV_sum
if nargin < 7, fit_propspread = 0; end
if nargin < 9, u_aversion = 0; end
if nargin < 10, saveresults = 1; end
if nargin < 11, graphics = 0; end

global rew_rng_state 
rew_rng_seed = 99;

results_dir = '/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out';
[~, str] = fileparts(data_file);
id = regexp(str,'(?<=fMRIEmoClock_)[\d_]+(?=_tc)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

%load data from CSV file
[data, y, u] = sceptic_get_data(data_file, n_steps);
[options, dim] = sceptic_get_options(data, nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps, u_aversion, graphics);

%Set up sigma noise for every point in u or hidden state?
rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
rew_rng_state=rng;
[~,idx] = unique(data.run);
conditions=data.rewFunc(idx);
sigma_noise = zeros(length(conditions), 1);
run_length = n_t/n_runs;
for i = 1:length(conditions)
    sni = repmat(std(arrayfun(@(x) RewFunction(x*100, conditions(i), 0), options.inF.tvec))^2, 1, run_length);
    sigma_noise(i) = mean(sni);
end
sigma_noise = mean(sigma_noise);
options.inF.sigma_noise = sigma_noise;

% Evolution function
h_name = @h_sceptic_fixed_decay_fmri;

% Observation function
g_name = @g_sceptic;

%this is not used
%onsets=NaN(1,size(y,2));
%for mm=1:size(y,2)
%    pos = find(y(:,mm)==1);
%    if isempty(pos), pos=NaN; end
%    onsets(mm) = pos;
%end

[posterior,out] = VBA_NLStateSpaceModel(y,u,h_name,g_name,dim,options);

%scepticrefit(posterior, out);

if saveresults
    %% save output figure
    % h = figure(1);
    % savefig(h,sprintf('results/%s_%s_multinomial%d_multisession%d_fixedParams%d', data_file,model,multinomial,multisession,fixed_params_across_runs))
    save(sprintf([results_dir, '/SHIFTED_U_CORRECT_%s_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d_sceptic_vba_fit'], id, model, multinomial,multisession,fixed_params_across_runs, u_aversion), 'posterior', 'out');
end
