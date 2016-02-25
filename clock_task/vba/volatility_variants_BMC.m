%% quick look at volatility variants

cd('low_res/')

modelnames = {'fixedUV'; 'sigVol'; 'sigVolLocal'; 'sigVolLocalNogamma'; 'sigVolLocPrecisionNogamma'};


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


for i = 1:4
    xlabel(out.options.handles.ha(i),'models')
    set(out.options.handles.ha(i),'xtick',1:length(modelnames), 'XTickLabel',char(modelnames))
end

h = figure(4);
file_str=input('What do you want save the figure as? ', 's');
saveas(h,[file_str])


%  without UV_sum

modelnames = {'sigVol'; 'sigVolLocal'; 'sigVolLocalNogamma'; 'sigVolLocPrecisionNogamma'};
L(1,:) = sigmavolatility.L;
L(2,:) = sig_vol_local.L;
L(3,:) = sig_vol_local_no_gamma.L;
L(4,:) = sig_vol_local_no_gamma_precision.L;

[posterior,out] = VBA_groupBMC(L);


for i = 1:4
    xlabel(out.options.handles.ha(i),'models')
    set(out.options.handles.ha(i),'xtick',1:length(modelnames), 'XTickLabel',char(modelnames))
end


