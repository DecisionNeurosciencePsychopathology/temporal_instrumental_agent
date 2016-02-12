%% quick look at volatility variants

cd('low_res/')

names = {'fixed_uv'; 'sigmavolatility'; 'sig_vol_local'; 'sig_vol_local_nogamma'; 'sig_vol_local_precision_nogamma'};


fixed_uv = load('grp_only_fixed_uv4_nbasis10_nsteps');
sigmavolatility = load('grp_only_kalman_sigmavolatility4_nbasis10_nsteps');
sig_vol_local = load('grp_only_kalman_sigmavolatility_local4_nbasis10_nsteps0_no_gamma');
sig_vol_local_no_gamma = load('grp_only_kalman_sigmavolatility_local4_nbasis10_nsteps1_no_gamma');
sig_vol_local_no_gamma_precision = load('grp_only_kalman_sigmavolatility_precision4_nbasis10_nsteps0_no_gamma');


L(1,:) = fixed_uv.L;
L(2,:) = sigmavolatility.L;
L(3,:) = sig_vol_local.L;
L(4,:) = sig_vol_local_no_gamma.L;
L(5,:) = sig_vol_local_no_gamma_precision.L;


options.families{1} = 1;
options.families{2} = 2:5;


[posterior,out] = VBA_groupBMC(L,options);
