
% calculates range of RTs to gauge the extent of exploration, both group and individual

close all
clear variables;

%% set up the group level

os = computer;


if strcmp(os(1:end-2),'PCWIN')
    behavfiles = glob('C:\kod\temporal_instrumental_agent\clock_task\subjects\*.csv');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'Alex')==1
        behavfiles = glob('/Users/localadmin/code/clock_smoothoperator/clock_task/subjects/*.csv');
        results_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
        group_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc';
    elseif strcmp(me(1:6),'dombax')==1
        behavfiles = glob('/Users/dombax/temporal_instrumental_agent/clock_task/subjects/*.csv');
        results_dir = '/Users/dombax/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
        group_dir = '/Users/dombax/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc';
    else
        behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
        addpath(genpath('/Users/dombax/temporal_instrumental_agent/clock_task/'));
        addpath(genpath('/Users/dombax/code/'));
    end
end

% get ID list

parfor sub = 1:length(behavfiles)
    str = behavfiles{sub};
    id = str2double(str(isstrprop(str,'digit')));
    fprintf('Fitting subject %d \r',sub)
    if strcmp(os(1:end-2),'PCWIN')
        data = readtable(sprintf('c:/kod/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
        vbadir = 'c:/kod/temporal_instrumental_agent/clock_task/vba';
        results_dir = 'E:\Users\wilsonj3\Google Drive\skinner\SCEPTIC\subject_fitting\vba_results\corrected_uv_sum';
    else
        [~, me] = system('whoami');
        me = strtrim(me);
        if strcmp(me,'Alex')==1
            data = readtable(sprintf('/Users/localadmin/code/clock_smoothoperator/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
            vbadir = '/Users/localadmin/code/clock_smoothoperator/clock_task/vba';
            results_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
            
        elseif strcmp(me(1:6),'dombax')==1
            data = readtable(sprintf('/Users/dombax/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
            vbadir = '/Volumes/bek/vba_results/uv_sum';
            results_dir = '/Users/dombax/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
            
        elseif strcmpi(me(1:14),'alexdombrovski')
            data = readtable(sprintf('/Users/alexdombrovski/code/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
            vbadir = '/Users/alexdombrovski/code/temporal_instrumental_agent/clock_task/vba';
            results_dir = '/Users/alexdombrovski/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
            
        else
            data = readtable(sprintf('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
        end
    end
    
    n_t = size(data,1);
    n_runs = n_t/50;
    trialsToFit = 1:n_t;
    
    
%     [~,idx] = unique(data.run);
%     conditions=data.rewFunc(idx);
%     
%     run_length = n_t/n_runs;
%     
    
    % Inputs
    rts(sub,:) = data{trialsToFit, 'rt'};
    
end



rts(rts==4000) = NaN;

rt_range.rts = rts;
rt_range.min_rts = min(rts');
rt_range.max_rts = max(rts');
rt_range.ranges = rt_range.max_rts - rt_range.min_rts;
rt_range.percent = rt_range.ranges./4000;
rt_range.percent_sd = var(rt_range.percent);
rt_range.percent_mean = mean(rt_range.percent);
rt_range.percent_min = min(rt_range.percent);
rt_range.percent_max = max(rt_range.percent);

n_runs = 8;

for ct = 1:n_runs;
block_rt_range.rts(:,:,ct) = rts(:,(ct-1)*50+1:ct*50);
block_rt_range.max_rts(:,ct) = max(block_rt_range.rts(:,:,ct),[],2);
block_rt_range.min_rts(:,ct) = min(block_rt_range.rts(:,:,ct),[],2);
block_rt_range.ranges(:,ct) = block_rt_range.max_rts(:,ct) - block_rt_range.min_rts(:,ct);
block_rt_range.percent(:,ct) = block_rt_range.ranges(:,ct)./4000;
block_rt_range.percent_sd(:,ct) = var(block_rt_range.percent(:,ct));
block_rt_range.percent_mean(:,ct) = mean(block_rt_range.percent(:,ct));
block_rt_range.percent_min(:,ct) = min(block_rt_range.percent(:,ct));
block_rt_range.percent_max(:,ct) = max(block_rt_range.percent(:,ct));
end

block_rt_range.subject_percent = mean(block_rt_range.percent,2);
block_rt_range.grp_mean_percent = mean(block_rt_range.subject_percent);
block_rt_range.grp_sd_percent = var(block_rt_range.subject_percent);

figure(1); clf;
boxplot(rt_range.percent);
ylabel('Range of RTs, proportion of interval');

h = figure(2); 
boxplot(block_rt_range.percent);
ylabel('Range of RTs, proportion of interval');
xlabel('Run, condition counterbalanced across subjects');
savefig(h,sprintf('proporion_of_interval_explored_by_subject_and_block'));

h = figure(3);
plot(block_rt_range.percent');

cd(group_dir);
%% save output
save(sprintf('rt_range'));
