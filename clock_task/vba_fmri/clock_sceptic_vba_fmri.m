function [posterior,out] = clock_sceptic_vba_fmri(data_file, so)
% fits SCEPTIC model to single-subject Clock Task subject data using VBA toolbox
% example call:
% [posterior,out]=clock_sceptic_vba_fmri('test_data_tcExport.csv', so)
% data_file:    name of CSV file containing raw trial-level data from clock task
% so: struct of options for fitting used in SCEPTIC model (setup by sceptic_validate_options)

global rew_rng_state 
rew_rng_seed = 99;

[~, str] = fileparts(data_file);
id = regexp(str,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)
so.id = id; %propagate into inF and inG for later extraction in group summaries

%load data from CSV file
[data, y, u] = sceptic_get_data(data_file, so);
[vba_options, dim] = sceptic_get_vba_options(data, so);

%Set up sigma noise for every point in u or hidden state?
rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
rew_rng_state=rng;
[~,idx] = unique(data.run);
conditions=data.rewFunc(idx);
sigma_noise = zeros(length(conditions), 1);
run_length = dim.n_t/vba_options.inF.n_runs;
for i = 1:length(conditions)
    sni = repmat(std(arrayfun(@(x) RewFunction(x*100, conditions(i), 0), vba_options.inF.tvec))^2, 1, run_length);
    sigma_noise(i) = mean(sni);
end
sigma_noise = mean(sigma_noise);
vba_options.inF.sigma_noise = sigma_noise;

[posterior,out] = VBA_NLStateSpaceModel(y, u, so.evo_fname, so.obs_fname, dim, vba_options);

posterior=add_transformed_params(posterior, so); %add transformed phi and mu into output objects

out.diagnostics = VBA_getDiagnostics(posterior, out); %pre-compute diagnostics in batch mode (gives Volterra outputs, param correlations, etc.)

if so.saveresults
  % save output figure
  % h = figure(1);
  % savefig(h,sprintf('results/%s_%s_multinomial%d_multisession%d_fixedParams%d', data_file,model,multinomial,multisession,fixed_params_across_runs))
  
  if ~isfield(so, 'output_dir')
    so.output_dir = '/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out';
  end
  
  save(sprintf([so.output_dir, '/sceptic_fit_%s_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d'], ...
    id, so.model, so.multinomial, so.multisession, so.fixed_params_across_runs, so.u_aversion), ...
    'posterior', 'out');
end
