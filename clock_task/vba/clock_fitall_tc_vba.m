allsubjs = dir('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/*.csv');
ids = NaN(length(allsubjs), 1);
parfor i = 1:length(allsubjs)
    ids(i) = str2double(regexprep(allsubjs(i).name,'fMRIEmoClock_(\d+)_tc_tcExport.csv','$1'));
    fprintf('processing id: %d\n', ids(i));
    %only run if no previous file
    if ~(exist(sprintf('results/%d_tc_vba_fit.mat', ids(i)), 'file') == 2)
        [posterior, out] = clock_tc_vba(ids(i), 0); %suppress figure
        [posterior, out] = clock_klambda_vba(ids(i), 0); %suppress figure
    else
        fprintf('exists: results/%d_tc_vba_fit.mat\n', ids(i));
    end
end