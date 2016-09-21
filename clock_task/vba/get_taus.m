% files = glob('E:/data/sceptic/vba_out/new_lambda_results/SHIFTED_U*kalman_uv*choice_tbf_autocorrel*.mat');
% 
% L = [];
% gamma = [];
% for i = 1:length(files)
%     load(files{i})
%     L(i) = out.F;
%     uv_sum_tau(i) = posterior.muTheta(1)./1000;
%     close all
% end
% 
% 
% files = glob('E:/data/sceptic/vba_out/new_lambda_results/SHIFTED_U*fixed_uv*choice_tbf_autocorrel*.mat');
% 
% L = [];
% gamma = [];
% for i = 1:length(files)
%     load(files{i})
%     L(i) = out.F;
%     fixed_uv_tau(i) = posterior.muTheta(1)./1000;
%     close all
% end



files = glob('E:/data/sceptic/vba_out/new_lambda_results/SHIFTED_U*fixed_decay_*_none_auto*.mat');
nbasis=24;
L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    S_decay(i,:)=out.suffStat.muX(nbasis+1,:);
    value_max_decay(i,:)=out.suffStat.muX(nbasis+2,:);
    rt_decay(i,:) = out.u(1,:);
    L_decay(i) = out.F;
end

files = glob('E:/data/sceptic/vba_out/new_lambda_results/SHIFTED_U*fixed_m*_none_auto*.mat');

L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    S_fixed(i,:)=out.suffStat.muX(nbasis+1,:);
    rt_fixed(i,:) = out.u(1,:);
    value_max_fixed(i,:)=out.suffStat.muX(nbasis+2,:);
    L_fixed(i) = out.F;
end