%loads in subjects' data and fits SCEPTIC models using VBA;
close all; clear;

%curpath = fileparts(mfilename('fullpath'));

% Note: this function looks for 'sceptic_dataset' and 'sceptic_model'
% as environment variables so that this script can be scaled easily for batch processing

so=[];

%so.model='fixed_uv_baked';
%so.model='decay_factorize_selective_psequate_fixedparams_fmri';
%so.model='decay_factorize_selective_psequate_fixedparams_meg';
so.mfx = 0; %so that output dir is set properly to have 'ffx' in path

%so.model='decay';
so = sceptic_validate_options(so); %initialize and validate sceptic fitting settings

%setup paths, parallel pools, etc. based on user environment
[so, poolobj, behavfiles] = sceptic_setup_environment(so);

%extract IDs for record keeping
[~,fnames]=cellfun(@fileparts, behavfiles, 'UniformOutput', false);
ids=cellfun(@(x) char(regexp(x,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match')), fnames, 'UniformOutput', false);

%puts a nested cell in each element
%ids=regexp(fnames,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

% Log evidence matrix
L = NaN(1,length(behavfiles));

% Subject statistics cell vector
s_all = cell(1,length(behavfiles));

% p = ProgressBar(length(behavfiles));

parfor sub = 1:length(behavfiles)
  fprintf('Fitting subject %d id: %s \n', sub, ids{sub});
  
  [posterior, out] = clock_sceptic_vba_fmri(behavfiles{sub}, so);
  s_all{sub} = extract_subject_statistics(posterior, out); %extract key statistics for each subject
  
  L(sub) = out.F;
  
  subj_id=ids{sub};

  %write out the basis as a csv file (for postprocessing) in one subject
  if sub == 1
    %save the basis as a csv file
    vnames = cellfun(@(x) strcat('Time_', num2str(x)), num2cell(1:so.ntimesteps), 'UniformOutput', false);
    basis_mat = array2table(out.options.inF.gaussmat, 'VariableNames', vnames);

    writetable(basis_mat, sprintf('%s/%s_%s_ffx_sceptic_basis.csv', so.output_dir, so.dataset, so.model));    
  end
  
  %parsave doesn't work in recent MATLAB versions...
  m=matfile(sprintf('%s/sceptic_fit_%s_%s_multinomial%d_multisession%d_fixedparams%d_uaversion%d', ...
    so.output_dir, ids{sub}, so.model, so.multinomial, so.multisession, ...
    so.fixed_params_across_runs, so.u_aversion), 'writable',true);
  
  m.posterior=posterior; m.out=out; m.subj_id=ids{sub}; m.subj_stats=s_all{sub};
  
  %parsave(sprintf('%s/sceptic_fit_%s_%s_multinomial%d_multisession%d_fixedparams%d_uaversion%d', ...
  %		  so.output_dir, ids{sub}, so.model, so.multinomial, so.multisession, ...
  %    so.fixed_params_across_runs, so.u_aversion), posterior, out);%, subj_id);
  
end

% p.stop;
delete(poolobj);

[group_global, group_trial_level] = extract_group_statistics(s_all, ...
  sprintf('%s/%s_%s_ffx_sceptic_global_statistics.csv', so.output_dir, so.dataset, so.model), ...
  sprintf('%s/%s_%s_ffx_sceptic_trial_outputs_by_timestep.csv', so.output_dir, so.dataset, so.model));

%save group outputs for now
save(sprintf('%s/group_fits_%s_%s', so.output_dir, so.model, so.dataset), 'ids', 'L', 'so', 's_all', 'group_global', 'group_trial_level');
