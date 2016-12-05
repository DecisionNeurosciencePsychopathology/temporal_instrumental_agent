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

os = computer;
if strcmp(os(1:end-2),'PCWIN')
    data_dir = 'E:/data/sceptic/vba_out/new_lambda_results/';
else
    data_dir = '/Volumes/bek/vba_results/';
end

multisession=1;

if multisession==0
    file_str_decay = 'SHIFTED_U*fixed_decay_*_none_auto*_total*.mat';
    file_str_fixed = 'SHIFTED_U*fixed_m*_none_auto*_total*.mat';
else
    file_str_decay = 'SHIFTED_U*fixed_decay_*_multisession1_*none_auto*_total*.mat';
    file_str_fixed = 'SHIFTED_U*fixed_m*_multisession1_*_none_auto*_total*.mat';
end

%Fixed decay
files = glob([data_dir file_str_decay]);
nbasis=24;
L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    value = out.suffStat.muX(1:nbasis,:);
    value = [zeros(nbasis,1) value(:,1:length(value)-1)];
    active_elements = value>0.01;
    gx = out.suffStat.gx;
    for j = 1:length(value)
        val_per_trial = value(:,j);
%         H_decay_shannon(i,j)=wentropy(val_per_trial(active_elements(:,j)),'shannon');
%         H_decay_gx_shannon(i,j) = wentropy(gx(:,j),'shannon');
%         H_decay_gx(i,j) = wentropy(gx(:,j),'logenergy');
%         
%         %Try normailizing it to be between 0 and 1?
%         H_decay_norm_gx_shannon(i,j) = wentropy((gx(:,j) - min(gx(:,j))) / ( max(gx(:,j)) - min(gx(:,j)) ),'shannon');
%         H_decay_norm_gx(i,j) = wentropy((gx(:,j) - min(gx(:,j))) / ( max(gx(:,j)) - min(gx(:,j)) ),'logenergy');
%         
%         
        
                
        %Lets try not using the wentropy
        my_shannon_decay(i,j) = calc_shannon_H(gx(:,j));
        
        %Let try to normailze the basis value basis funcitons into pdf
        val_p = val_per_trial/sum(val_per_trial); %so it sums to 1
        my_shannon_decay_val(i,j) = calc_shannon_H(val_p);
        wentropy_shannon_decay_val(i,j)  = wentropy(val_p,'shannon');
        
    end
    
    H_decay(i,:)=out.suffStat.muX(nbasis+1,:);
    value_max_decay(i,:)=out.suffStat.muX(nbasis+2,:);
    rt_decay(i,:) = out.u(1,:);
    L_decay(i) = out.F;
end


%Fixed
files = glob([data_dir file_str_fixed]);

L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    
    value = out.suffStat.muX(1:nbasis,:);
    value = [zeros(nbasis,1) value(:,1:length(value)-1)];
    active_elements = value>0.01;
    gx = out.suffStat.gx;
    for j = 1:length(value)
        val_per_trial = value(:,j);
%         H_fixed_shannon(i,j)=wentropy(val_per_trial(active_elements(:,j)),'shannon');
%         H_fixed_gx_shannon(i,j) = wentropy(gx(:,j),'shannon');
%         H_fixed_gx(i,j) = wentropy(gx(:,j),'logenergy');
%         
%         
%         %Try normailizing it to be between 0 and 1?
%         H_fixed_norm_gx_shannon(i,j) = wentropy((gx(:,j) - min(gx(:,j))) / ( max(gx(:,j)) - min(gx(:,j)) ),'shannon');
%         H_fixed_norm_gx(i,j) = wentropy((gx(:,j) - min(gx(:,j))) / ( max(gx(:,j)) - min(gx(:,j)) ),'logenergy');
%         
        %Lets try not using the wentropy
        my_shannon_fixed(i,j) = calc_shannon_H(gx(:,j));
        
        
        %Let try to normailze the basis value basis funcitons into pdf
        val_p = val_per_trial/sum(val_per_trial); %so it sums to 1
        my_shannon_fixed_val(i,j) = calc_shannon_H(val_p);
        wentropy_shannon_fixed_val(i,j)  = wentropy(val_p,'shannon');
        
        
    end
    H_fixed(i,:)=out.suffStat.muX(nbasis+1,:);
    rt_fixed(i,:) = out.u(1,:);
    value_max_fixed(i,:)=out.suffStat.muX(nbasis+2,:);
    L_fixed(i) = out.F;
end