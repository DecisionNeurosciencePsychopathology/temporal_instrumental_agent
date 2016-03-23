function subj_data=analyze_subjects_rt_behavior
%Find global max and dislodgment statistics for each subject

behavfiles = glob('../subjects/*.csv');
graphics=0;


% get ID list
id = NaN(length(behavfiles),1);

%Run through all subjects
for sub = 1:length(behavfiles)
    str = behavfiles{sub};
    id(sub) = str2double(str(isstrprop(str,'digit')));
    fprintf('Analyzing subject %d id: %d \r',sub, id(sub))
    subj_data.(['id' num2str(id(sub))]) =...
        subject_rts_distribution_sub_proc(id(sub),graphics);
end

%Might need to change up the data structure a bit but this works for now I
%think...

Snames = fieldnames(subj_data);
for loopIndex = 1:numel(Snames)
    global_max_hit(loopIndex,1)=...
        subj_data.(Snames{loopIndex}).global_max_reached;
    dislodge_hit(loopIndex,1)=...
        subj_data.(Snames{loopIndex}).dislodgment_occurrence;
    no_recover(loopIndex,1)=...
        subj_data.(Snames{loopIndex}).no_recover;
end

global_max_ratio = mean(global_max_hit~=0);
dislodgment_ratio = mean(dislodge_hit~=0);
per_run_global_max_ratio = (global_max_hit/6);
per_run_dislodgment_ratio = (dislodge_hit/6);
no_recover_ratio = mean(no_recover);
fprintf('\n\nProportion of global max achieved per subject: %.2f%%\n',global_max_ratio*100);
fprintf('Proportion of dislodgment per subject: %.2f%%\n',dislodgment_ratio*100);

fprintf('Proportion of global max achieved per subject per run: %.2f runs\n',mean(global_max_hit));
fprintf('Proportion of dislodgment per subject per run: %.2f runs\n',mean(dislodge_hit));

fprintf('Proportion of subjects who did not recover after dislodging: %.2f%%\n',no_recover_ratio*100)

figure(8)
histogram(per_run_global_max_ratio)
figure(80)
histogram(per_run_dislodgment_ratio)
