function subj_data=subject_rts_distribution_sub_proc
%Determining the global max distributions for each subject
%
id = 10637; %make this an input later
close all;
%% Grab data
data = readtable(sprintf('/Users/dombax/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);

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
graphics = 1; %Make an arg
reach_criteria = 6;
fprintf('Subject: %d\n',id);

for i = 1:length(subj_rts)
    if isempty(subj_rts{i}), continue; end %Skip non IEV|DEV runs
    
    %What condition is it
    current_cond = unique(data.rewFunc(data.run==i));
    
    %Initalize struct for data storage
    %subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)])=[]; 
    
    %Determine if subejct meets criteria for this run (i.e has 4 RTs in
    %target window)
    if strcmp('IEV',current_cond)
        target_window = 3000:4000; %Data should always be in ms
        subj_data.(['id' num2str(id)]).rtsIEV(i,:)=...
            subj_rts{i};
    else
        target_window = 250:1250; %Should this be 0 -> 1000?
        subj_data.(['id' num2str(id)]).rtsDEV(i,:)=...
            subj_rts{i};
    end
    
    criteria_idx=ismember(subj_rts{i},target_window)';
    [target_event, trials_in_target_zone] = event_locator(criteria_idx);
    max_target_run=max(target_event);
    
    %Find dislodgement points
    dislodgement_idx = false(1,length(subj_rts{1})); %Might just grab num trials instead?
    for k = reach_criteria+1:length(criteria_idx)
        if (~criteria_idx(k)) && (dislodgement_idx(k-1) || sum(criteria_idx(k-reach_criteria:k-1))==reach_criteria)
            dislodgement_idx(k) = true;
        end
    end
            
    [dislodge_event, dislogde_trial_count] = event_locator(dislodgement_idx); %Total number of trials that were considered dislodged
    max_dislodgement_run=max(dislodge_event); %Highest dislodge run
    dislodge_occurrence = find( (dislodgement_idx > 0) );
    dislodge_occurrence =dislodge_occurrence(find(diff([0 dislodge_occurrence])-1));  %#ok<FNDSB> %Where the dislodgements occured
    total_num_dislodgments = numel(dislodge_occurrence); %Total number of times subject dislodged from target 
    %cellfun(@ismember,target_window)
    if ge(max_target_run,reach_criteria)
        fprintf('Global max hit for run %d, %s\n',i,current_cond{:})
    end
    
    %Make some pretty plots of RTs
    x = 1:50;
    if graphics
        figure(i+10)
        clf;
        subplot(2,1,1)
        scatter(1:length(subj_rts{i}),subj_rts{i});
        hold on;
        scatter(x(criteria_idx),subj_rts{i}(criteria_idx),'filled');
        axis([0 50 250 4000])
        title(['Run: ', num2str(i), 'Condition: ', current_cond])
        xlabel('Trial');
        ylabel('Time');
        subplot(2,1,2)
        edges = [0 1001 2001 3001 4000];
        histogram(subj_rts{i},edges);
        xlabel('Time');
        ylabel('Bin');
    end    
    h = figure(i+10);
    
    %Other stats
    out_of_range_rts = subj_rts{i}(~criteria_idx);
    dislodgement_rts = subj_rts{i}(dislodgement_idx);
    %Determine when subject was out of target window and compute the
    %average
    mean_out_of_target_range_rt = mean(abs(out_of_range_rts-max(target_window))); %PERHAPS rewrite this 'subj_rts{i}(~criteria_idx)' as out of bounds rts?
    max_dislodgemnet_time = max(dislodgement_rts); %Max dislodgement run timepoint
    mean_dislodgement_rt =  mean(abs(dislodgement_rts)-max(target_window));
    
    
    
    %Put stat in a struct although maybe we should put everything into a
    %matrix and create fieldnames and use setfield? I don't like this kind
    %of syntax...
    subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)]).max_target_run=max_target_run;
    subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)]).trials_in_target_zone=trials_in_target_zone;
    subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)]).max_dislodgement_run=max_dislodgement_run;
    subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)]).total_num_dislodgments=total_num_dislodgments;
    subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)]).dislogde_trial_count=dislogde_trial_count;
    subj_data.(['id' num2str(id)]).([current_cond{:} num2str(i)]).figure_data=h;
    %Add total off target, max dislodgement difference, maybe include a raw
    %which contains the indexs of dislodge and target.
    

end
 %Perhaps make a wait for button if id is on a certain list or meets a certain criteria to investigate further?
 
 %Get the rts the subject chose for each run
 subj_data.(['id' num2str(id)]).rtsIEV =...
     remove_zeros(subj_data.(['id' num2str(id)]).rtsIEV);
  subj_data.(['id' num2str(id)]).rtsDEV =...
     remove_zeros(subj_data.(['id' num2str(id)]).rtsDEV);
 
    function [events, count]=event_locator(idx)
        events=diff( [0 (find( ~ (idx > 0) ) )...
        numel(idx) + 1] ) - 1;
        
        count = sum(idx);
        
        function s=remove_zeros(s)
            s(s==0)=[];
        
 