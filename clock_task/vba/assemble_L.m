% The script Jon and Alex used for final group comparison on 02/11/16, G-d
% help them

L8 = load('L8_fixed_ksoftmax_kprocnoise_kuvsum_ksigvol_kuvlog_fixeduv_q');
Lk = L8.L8(2:6,:);
Lf = L8.L8([1 7],:);
load('grp_qlearning_nsteps_40_stepWise.mat');
Lq = L;
Lall = [Lf;Lk;Lq];
options.families{1} = [1 2];
options.families{2} = [3 4 5 6 7];
options.families{3} = [8];
modelnames = {'fixed', 'fixedUv', 'kSoftmax', 'kProcNoise', 'kUvSum', 'kSigVol', 'kUvLog', 'Q'};
[posterior,out] = VBA_groupBMC(Lall,options);

for i = [1:length(out.options.handles.ha)-2 7]
xlabel(out.options.handles.ha(i),'models')
set(out.options.handles.ha(i),'xtick',1:length(modelnames), 'XTickLabel',char(modelnames))
end
set(out.options.handles.ha(6),'xtick',1:3, 'XTickLabel',char({'fixed', 'kalman', 'q-learning'}));

%% CHECK!
fignum = 10;
h = figure(fignum);
savefig(h,'grp_BMC_final_USE_THIS_ONE');