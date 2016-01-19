allsubjs = dir('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/*.csv');
ids = NaN(length(allsubjs), 1);
parfor i = 1:length(allsubjs)
    ids(i) = str2double(regexprep(allsubjs(i).name,'fMRIEmoClock_(\d+)_tc_tcExport.csv','$1'));
    [posterior, out] = clock_tc_vba(ids(i), 0); %suppress figure
end