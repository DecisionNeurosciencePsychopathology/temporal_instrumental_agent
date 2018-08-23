function [so] = sceptic_validate_options(so)
  % nbasis:      8 works well, 4 also OK
  % multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
  % multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
  % fixed_params_across_runs -- self-explanatory
  % fit_propspread -- makes temporal generalization within the eligibility trace a free parameter
  % ntimesteps:      number of time bins
  % u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum

  if nargin < 1, so=[]; end

  %if user specifies a dataset upstream, don't check environment variable
  if isfield(so, 'sceptic_dataset')
    fprintf('Using user-specified dataset: %s.\n   Will not check sceptic_dataset environment variable', so.sceptic_dataset);
  else
    so.dataset=getenv('sceptic_dataset');
    if strcmpi(so.dataset, '')
      so.dataset='mmclock_fmri'; %or mmclock_fmri or specc_fmri
    end
  end
  
  %same principle for model variant
  if isfield(so, 'sceptic_model')
    fprintf('Using user-specified model: %s.\n   Will not check sceptic_model environment variable', so.sceptic_model);
  else
    so.model=getenv('sceptic_model');
    if strcmpi(so.model, '')
      so.model = 'decay'; % will run to get value and prediction errors.
    end
  end
  
  if ~isfield(so, 'nbasis'), so.nbasis=24; end
  if ~isfield(so, 'multinomial'), so.multinomial=1; end
  if ~isfield(so, 'multisession'), so.multisession=0; end
  if ~isfield(so, 'fixed_params_across_runs'), so.fixed_params_across_runs=1; end
  if ~isfield(so, 'fit_propspread'), so.fit_propspread=0; end
  if ~isfield(so, 'max_prop_spread'), so.max_prop_spread = .0125; end %used for fixed propspread case
  if ~isfield(so, 'ntimesteps'), so.ntimesteps=40; end
  if ~isfield(so, 'u_aversion'), so.u_aversion=1; end  % default based on Cognition analyses
  if ~isfield(so, 'graphics'), so.graphics=0; end
  if ~isfield(so, 'model'), so.model='fixed'; end
  if ~isfield(so, 'trials_per_run'), so.trials_per_run=50; end
  if ~isfield(so, 'trial_length'), so.trial_length=4000; end %trial length in ms (used for discretizing to bins)
  if ~isfield(so, 'saveresults'), so.saveresults=-0; end %used in some places to denote whether to save fitted outputs

  if strcmpi(so.model,'fixed')
    so.evo_fname = @h_sceptic_fixed_fmri;
    so.hidden_states=2; %track basis weights (value), as well as PE as tag-along state
    so.n_theta=1; %learning rate
  elseif strcmpi(so.model,'decay')
    so.evo_fname = @h_sceptic_fixed_decay_fmri;
    so.hidden_states=3; %track basis weights (value), as well as PE and decay as tag-along states
    so.n_theta=2; %learning rate and selective maintenance parameters (decay outside of the eligibility trace)
  elseif strcmpi(so.model,'fixed_UV')
    so.evo_fname = @h_sceptic_kalman_fmri;
    so.hidden_states=3; %track basis weights (value), as well as PE and uncertainty
    so.n_theta=2; %learning rate alpha and uncertainty sensitivity parameter tau

  end
  
  so.obs_fname = @g_sceptic; %just regular softmax observation rule for now
  
  if strcmpi(so.dataset, 'mmclock_meg')
    so.trials_per_run=63; %63 for MEG, 50 for fMRI
  else
    so.trials_per_run=50;
  end
  
end
