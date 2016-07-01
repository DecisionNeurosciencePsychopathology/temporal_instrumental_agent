
modelnames = {'fixed' 'fixed_uv' 'fixed_decay' 'kalman_softmax' 'kalman_processnoise' 'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'};
% results_dir, '/SHIFTED_U_CORRECT%d_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d_sceptic_vba_fit_fixed_prop_spread'], id, model, multinomial,multisession,fixed_params_across_runs, u_aversion)
% files = glob('/Volumes/bek/vba_results/sigma_volatility_variants/precision/*.mat');

for j = 1:length(modelnames)
    
    file_str = ['/Volumes/bek/vba_results/*_', modelnames{j}, '_m*_fixed_prop_spread.mat']; 
     files = glob(file_str);
    for i = 1:length(files)
        load(files{i})
        L(j,i) = out.F;
    end
end