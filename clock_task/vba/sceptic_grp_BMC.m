%runs group model comparison on SCEPTIC log model evidence data from sceptic_fit_group_vba.m
%Your data file should contain the subjects' posterior log evidences as
%well the names of the model ran


close all;

%Set path to
f_name = ''; %Set full path to posterior log evidences (We use the variable L)

%Load in proper data file
group_data=load([f_name]);

%Do you want to run models in families?
family_flag =1;

%% really manual 2-model comparison
sceptic.L = group_data.L; %Rename L as appropriate
modelnames = group_data.modelnames; %Rename modelnames as appropriate

if family_flag
    % Define 'fixed' as the first family, 'kalman', as the second, q third
    idx = strfind(modelnames, 'fixed');
    options.families{1} = find(not(cellfun('isempty', idx))); 
    idx = strfind(modelnames, 'kalman');
    options.families{2} = find(not(cellfun('isempty', idx)));
    [posterior,out] = VBA_groupBMC(sceptic.L,options);
else
    [posterior,out] = VBA_groupBMC(sceptic.L);
end


%% Clean up model strings to make graphic look pretty
% modelnames = cellfun(@strrep, modelnames', repmat({'_'},length(modelnames),1), repmat({' '},length(modelnames),1), 'UniformOutput', false);
for i = 1:length(out.options.handles.ha)-2
    xlabel(out.options.handles.ha(i),'models')
    set(out.options.handles.ha(i),'xtick',1:length(modelnames), 'XTickLabel',char(modelnames))
end


%% save output figure
h = figure(1);
file_str=input('What do you want save the figure as? ', 's');
saveas(h,[file_path file_str])