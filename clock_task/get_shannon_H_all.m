%% Estimate Shannon entropy for all subjects


%Discalimer!! This code is terrible are needs rewritten if it will ever be
%shared!

multisession=0;

os = computer;
if strcmp(os(1:end-2),'PCWIN')
    data_dir = 'E:/data/sceptic/vba_out/new_lambda_results/';
else
    data_dir = '/Volumes/bek/vba_results/';
end

if multisession==0
    %file_str_decay = 'SHIFTED_U*fixed_decay_*_none_auto*_total*.mat';
    %file_str_fixed = 'SHIFTED_U*fixed_m*_none_auto*_total*.mat';
    file_str_decay = '*meg*.mat';
    file_str_fixed = 'SHIFTED_U*fixed_m*.mat';
    %This was added later when we needed to change the entropy function to
    %utilize the hist f(x)
    %data_dir = '/Volumes/bek/vba_results/unisession_entropy_mar_6_2017/';
    %data_dir = '/Volumes/bek/vba_results/rand_priors_unisession/';
    data_dir = 'E:/data/sceptic/meg_clock/individual_results/fixed_decay/';
else
    %     file_str_decay = 'SHIFTED_U*fixed_decay_*_multisession1_*none_auto*_total*.mat';
    %     file_str_fixed = 'SHIFTED_U*fixed_m*_multisession1_*_none_auto*_total*.mat';
    
    
    %file_str_decay = 'SHIFTED_U*fixed_decay_*_multisession1_*.mat';
    file_str_decay = '*meg*.mat';
    file_str_fixed = 'SHIFTED_U*fixed_m*_multisession1_*.mat';
    
    %This was added later when we needed to change the entropy function to
    %utilize the hist f(x)
    %data_dir = '/Volumes/bek/vba_results/multisession_1_mar_7_2017_all_models/';
    data_dir = 'E:/data/sceptic/meg_clock/individual_results/fixed_decay/';
end

%Fixed decay
files = glob([data_dir file_str_decay]);
%These could be pulled from out also
nbasis=24;
nruns = 8;
L = [];
gamma = [];
%% PDF that goes into the entropy calculation: TBF element-wise value?
for i = 1:length(files)
    %Because I don't want to rewrite the entire script yet, as I'm not sure
    %if we will use it much after making plots, just save the data you want
    %to load in, then reload it for subseqent executions of this script.
    if (~exist('decay_data','var')) || length(decay_data)~=length(files)
        decay_data(:,:,i)=load(files{i});
    end
    
    %Load in data
    out = decay_data(i).out;
 
    %Need to reshape the value for multisession
    if multisession
        %We may have more hidden state than just the basis funcitons, in
        %order to calculate this just keep dividing the rows of muX by
        %nbasis + 1 until you get nruns (which should be 8)
        ct=0;
        [r,~]=size(out.suffStat.muX);
        while  r/(nbasis+ct)~=nruns
            ct = ct+1;
        end
        mux_idx = nbasis + ct;
        value = [];
        i1=1;
        i2=nbasis;
        run_len = 1:50;
        for ii = 1:nruns
            tmp_val = out.suffStat.muX(i1:i2,run_len);
            value = [value tmp_val];
            i1 = i1 + mux_idx;
            i2 = i2 + mux_idx;
            run_len = run_len +50;
        end
    else
        value = out.suffStat.muX(1:nbasis,:);
    end
    
    %Just shifting by 1 col (why again?)
    value = [zeros(nbasis,1) value(:,1:length(value)-1)];
    %     active_elements = value>0.01;
    %     gx = out.suffStat.gx;
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
        %         my_shannon_decay(i,j) = calc_shannon_H_hist(gx(:,j));
        %my_shannon_decay(i,j) = calc_entropy_modd(gx(:,j)'); %Transpose only for modd version!!!
        
        %Normailze the TBF value funciton into pdf
        val_p = val_per_trial/sum(val_per_trial); %so it sums to 1
        my_shannon_decay_val(i,j) = calc_shannon_H(val_p);
        %my_shannon_decay_val(i,j) = calc_entropy_modd(val_p');
        %wentropy_shannon_decay_val(i,j)  = wentropy(val_p,'shannon');
        
    end
    
    %H_decay(i,:)=out.suffStat.muX(nbasis+1,:);
    %value_max_decay(i,:)=out.suffStat.muX(nbasis*2+2,:);
    %rt_decay(i,:) = out.u(1,:);
    %L_decay(i) = out.F;
end


%Fixed
files = glob([data_dir file_str_fixed]);

L = [];
gamma = [];
for i = 1:length(files)
    if (~exist('fixed_data','var')) || length(fixed_data)~=length(files)
        fixed_data(:,:,i)=load(files{i});
    end
    
    %Load in data
    out = fixed_data(i).out;
 
    %Need to reshape the value for multisession
    if multisession
        %We may have more hidden state than just the basis funcitons, in
        %order to calculate this just keep dividing the rows of muX by
        %nbasis + 1 until you get nruns (which should be 8)
        ct=0;
        [r,~]=size(out.suffStat.muX);
        while  r/(nbasis+ct)~=nruns
            ct = ct+1;
        end
        mux_idx = nbasis + ct;
        value = [];
        i1=1;
        i2=nbasis;
        run_len = 1:50;
        for ii = 1:nruns
            tmp_val = out.suffStat.muX(i1:i2,run_len);
            value = [value tmp_val];
            i1 = i1 + mux_idx;
            i2 = i2 + mux_idx;
            run_len = run_len +50;
        end
    else
        value = out.suffStat.muX(1:nbasis,:);
    end
    
    
    value = [zeros(nbasis,1) value(:,1:length(value)-1)];
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
        %         my_shannon_fixed(i,j) = calc_shannon_H(gx(:,j));
        
        
        % Normailze the TBF value funciton into pdf
        val_p = val_per_trial/sum(val_per_trial); %so it sums to 1
        %         my_shannon_fixed_val(i,j) = calc_shannon_H(val_p);
        my_shannon_fixed_val(i,j) = calc_shannon_H(val_p);
        
        %my_shannon_fixed_val(i,j) = calc_entropy_modd(val_p');
        %wentropy_shannon_fixed_val(i,j)  = wentropy(val_p,'shannon');
        
        
    end
    %H_fixed(i,:)=out.suffStat.muX(nbasis+1,:);
    rt_fixed(i,:) = out.u(1,:);
    value_max_fixed(i,:)=out.suffStat.muX(nbasis+2,:);
    L_fixed(i) = out.F;
end