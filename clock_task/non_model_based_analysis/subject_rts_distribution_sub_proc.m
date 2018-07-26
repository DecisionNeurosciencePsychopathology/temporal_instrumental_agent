function subj_data=subject_rts_distribution_sub_proc(id,graphics)
%Determining the global max distributions for each subject
%
%id = 10637; %make this an input later
close all;
%% Grab data
data = readtable(sprintf('~/code/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);

%% Determine which runs are IEV or DEV and grab the RTs
runs=unique(data.run);
run_of_interest = zeros(length(runs),1);
subj_rts = cell(length(runs),1);
for i = 1:length(runs)
    run_of_interest(i,1) = ismember('IEV',unique(data.rewFunc(data.run==runs(i))))...
        | ismember('DEV',unique(data.rewFunc(data.run==runs(i)))) ;
    
    if logical(run_of_interest(i,1))
        temp_rts = data.rt(data.run==i);
        %subj_rts{i,1} = temp_rts(end-19:end); %only use last 20 trials
        subj_rts{i,1} = temp_rts;
    end
end

%% Evaluate the distribution of RTs
% Max for IEV is max of ntimesteps. Max for DEV is 1.5 sec
reach_criterion = 6; %How many decisions subject needs to make in a row
n=50; %number of trials
%Initialize variables
global_run_ct=0;
dislodge_run_ct=0;
no_recover_ct = 0;
subj_data.global_max_reached=0;
subj_data.dislodgment_occurrence=0;
subj_data.no_recover = 0;

%Start main for loop
for i = 1:length(subj_rts)
    if isempty(subj_rts{i}), continue; end %Skip non IEV|DEV runs
    
    %What condition is it
    current_cond = unique(data.rewFunc(data.run==i));
    
    %Determine if subejct meets criteria for this run (i.e has 4 RTs in
    %target window)
    if strcmp('IEV',current_cond)
        target_window = 3500:4000; %Data should always be in ms
        subj_data.rtsIEV(i,:)=...
            subj_rts{i};
    else
        target_window = 250:750; %Should this be 0 -> 1000?
        subj_data.rtsDEV(i,:)=...
            subj_rts{i};
    end
    
    
    criteria_idx=ismember(subj_rts{i},target_window)'; %What trials were in the target window
    target_data = event_locator(criteria_idx,reach_criterion); %Compute the target zone data
    max_target_run=max(target_data.max);     %The max run in the window
    total_target_runs_this_run = sum(target_data.total); %Total number of times subject went on a run in the target zone
    
    %Find dislodgment points
    dislodgment_idx = false(1,n); %Might just grab num trials instead?
    for k = reach_criterion+1:length(criteria_idx)
        if (~criteria_idx(k)) && (dislodgment_idx(k-1) || sum(criteria_idx(k-reach_criterion:k-1))==reach_criterion)
            dislodgment_idx(k) = true;
        end
    end
    
    
    dislodge_data = event_locator(dislodgment_idx,reach_criterion); %Generate dislodge data
    max_dislodgment_run=dislodge_data.max; %Highest dislodge run
    total_num_dislodgments_this_run = sum(dislodge_data.total); %Total number of times subject dislodged from target this run
    
    
    %Counting global max reached over runs
    if ge(max_target_run,reach_criterion)
        %fprintf('Global max hit for run %d, %s\n',i,current_cond{:})
        global_run_ct = global_run_ct+1;
        subj_data.global_max_reached = ...
            global_run_ct;
    end
    
    %Counting dislodgments overs runs
    if ge(max_dislodgment_run,reach_criterion)
        %fprintf('dislodgment occurred for run %d, %s\n',i,current_cond{:})
        dislodge_run_ct = dislodge_run_ct + 1;
        subj_data.dislodgment_occurrence = ...
            dislodge_run_ct;
        
        %Additionally does the subject ever recover after last
        %dislodgment?
        if max(dislodge_data.idx_and_length(:,1))>max(target_data.idx_and_length(:,1))
            no_recover_ct = no_recover_ct + 1;
            subj_data.no_recover = ...
            no_recover_ct;
        end
            
    end
    
    %Count climbs
    %     b=diff(subj_rts{i});
    %     rt_decreases = b<0;
    %     climbs = find(conv(double(rt_decreases==0),ones(1,reach_criterion),'valid')==reach_criterion); %// find k zeros
    %     if isempty(climbs)
    %         climb_count(i) = 0;
    %     else
    %     climbs = climbs([true diff(climbs)>n]); %// remove similar indices, indicating n+1, n+2... zeros
    %     climb_count(i) = length(climbs);
    %     climb_trials(i,:) = climbs:climbs+reach_criterion;
    %     end
    
    
    %Make some pretty plots of RTs
    if logical(graphics)
        h=make_a_plot(subj_rts{i},current_cond,i,criteria_idx);
        %Save figure
        subj_data.([current_cond{:} num2str(i)]).figure_data=h;
    end
    
    
    %Put stat in a struct although maybe we should put everything into a
    %matrix and create fieldnames and use setfield? I don't like this kind
    %of syntax...
    subj_data.([current_cond{:} num2str(i)]).max_target_run=max_target_run;
    %subj_data.([current_cond{:} num2str(i)]).trials_in_target_zone=trials_in_target_zone;
    subj_data.([current_cond{:} num2str(i)]).max_dislodgment_run=max_dislodgment_run;
    subj_data.([current_cond{:} num2str(i)]).total_num_dislodgments_this_run=total_num_dislodgments_this_run;
    subj_data.([current_cond{:} num2str(i)]).total_target_runs_this_run = total_target_runs_this_run;
    %subj_data.([current_cond{:} num2str(i)]).dislogde_trial_count=dislogde_trial_count;
    %Add total off target, max dislodgment difference, maybe include a raw
    %which contains the indexs of dislodge and target.
    
end

%Get the rts the subject chose for each run
subj_data.rtsIEV =...
    remove_zeros(subj_data.rtsIEV);
subj_data.rtsDEV =...
    remove_zeros(subj_data.rtsDEV);

function out=event_locator(idx,reach_crit)
%Clean up this funciton..
%Find where the events occured
events=diff( [0 (find( ~ (idx > 0) ) )...
    numel(idx) + 1] ) - 1;

%Get a total count
count = sum(idx);

%Find and index the events
event_occurrence = find( (idx > 0) );
event_occurrence_idx =event_occurrence(find(diff([0 event_occurrence])-1));
a=ismember(event_occurrence,event_occurrence_idx);
event_length=diff([(find(a)) numel(a) + 1] );
out.events = events;
out.count=count;
out.idx_and_length = [event_occurrence_idx; event_length]';
if size(out.idx_and_length,2)>=2
    out.max = max(out.idx_and_length(:,2));
    out.total = sum(out.idx_and_length(:,2)>=reach_crit);
else
    out.max = 0;
    out.total = 0;
end

function s=remove_zeros(s)
s(s==0)=[];

function h=make_a_plot(rts,curr_cond,idx,crit_idx)
h=figure(idx+10);
clf;
x=1:50;
subplot(2,1,1)
scatter(1:length(rts),rts);
hold on;
scatter(x(crit_idx),rts(crit_idx),'filled');
axis([0 50 250 4000])
title(['Run: ', num2str(idx), 'Condition: ', curr_cond])
xlabel('Trial');
ylabel('Time');
subplot(2,1,2)
edges = [0 1001 2001 3001 4000];
histogram(rts,edges);
xlabel('Time');
ylabel('Bin');


%% Code we may come back to
%Other stats
%     out_of_range_rts = subj_rts{i}(~criteria_idx);
%     dislodgment_rts = subj_rts{i}(dislodgment_idx);
%Determine when subject was out of target window and compute the
%average
%     mean_out_of_target_range_rt = mean(abs(out_of_range_rts-max(target_window))); %PERHAPS rewrite this 'subj_rts{i}(~criteria_idx)' as out of bounds rts?
%     max_dislodgemnet_time = max(dislodgment_rts); %Max dislodgment run timepoint
%     mean_dislodgment_rt =  mean(abs(dislodgment_rts)-max(target_window));

%     dislodge_occurrence = find( (dislodgment_idx > 0) );
%     dislodge_occurrence =dislodge_occurrence(find(diff([0 dislodge_occurrence])-1));  %#ok<FNDSB> %Where the dislodgments occured
%     total_num_dislodgments = numel(dislodge_occurrence); %Total number of times subject dislodged from target
