%Script to optimize parameters of skeptic model for various subjects  over
%runs, i.e. we are optimizing over runs, not over subjects.

global behavfiles sub model_names i


%load in behav files
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    behavfiles = glob('C:\kod\temporal_instrumental_agent\clock_task\subjects\*.csv');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'Alex')==1
        behavfiles = glob('/Users/localadmin/clock_smoothoperator/clock_task/subjects/*.csv');
    else
        behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
    end
end



%To start let's use table to create the data hub for softmax, v-discounted,
%and uv-discounted, to store our variables, initials values, and bounds.
model_names = {'value_softmax'; 'v_discounted'; 'uv_discounted'};
%Initial paramters
%1-value_softmax prop_spread, beta (temperature)
%2-v_discounted prop_spread, beta (temperature), kappa (discount factor), lambda (Go NoGo influence factor)
%3-uv_discounted prop_spread, beta (temperature), kappa (discount factor), tau (V vs. U weight)
if ~exist('default_info.mat','file')
%     initial_params = { [.2261 .01]; [.2261 .01 .005 .01]; [.2261 .01 .005 .95] };
%     lower_bounds = { [0.01 .001]; [0.01 .001 -.001 -.001]; [0.01 .001 -.001 .9] };
%     upper_bounds = { [0.06 1]; [0.06 1 .01 .01]; [0.06 1 .01 1] };
%JW: 9/10/15 changed betas, the small range was causing infs to appear
    initial_params = { [1]; [1 .005 .005]; [1 .005 .95] };
    lower_bounds = { [.01]; [.01 -.1 -.1]; [.01 -1 .5] };
    upper_bounds = { [10000]; [10000 .1 .1]; [10000 1 1] };
    %num_vars = [2; 4; 4];
    num_vars = [1; 3; 3];
    default_info= table(initial_params, lower_bounds, upper_bounds,...
        num_vars,'RowNames', model_names);
    save default_info default_info
else
    load('default_info.mat');
end



%Create main data table if it doesn't exist else just load it
if ~exist('subj_opt_vars.mat','file') && exist('default_info','var')
    idchars = regexp(behavfiles,'\d','match');
    ids=cellfun(@(x)(cell2mat(x)),idchars, 'UniformOutput',false);
    value_softmax = zeros(length(ids),length(default_info.initial_params{1}));
    value_softmaxFga = zeros(length(ids),1);
    v_discounted = zeros(length(ids),length(default_info.initial_params{2}));
    v_discountedFga = zeros(length(ids),1);
    uv_discounted = zeros(length(ids),length(default_info.initial_params{3}));
    uv_discountedFga = zeros(length(ids),1);
    subj_opt_vars=table(value_softmax, value_softmaxFga, v_discounted, v_discountedFga,...
        uv_discounted, uv_discountedFga,'RowNames',ids);
    save subj_opt_vars subj_opt_vars    
    clear value_softmax v_discounted uv_discounted; %Remove vars if needed
else
    load('subj_opt_vars.mat')
end

% main shell loop
for i = 2:numel(model_names) %model loop...9/10/15 since prop spread is the major player in the softmax version let's come back to it.
    for sub = 1:height(subj_opt_vars) %subject loop
        [subj_opt_vars.(model_names{i})(sub,:), subj_opt_vars.([model_names{i} 'Fga'])(sub)] = GA_skeptic_optimize(default_info.initial_params{i},...
            default_info.lower_bounds{i},default_info.upper_bounds{i},default_info.num_vars(i));
        
        %Everytime the model completes for one subject lets save the table and export the vars
        %to a text document for real time analysis
        save subj_opt_vars subj_opt_vars
        for foo=1:height(subj_opt_vars)
            writetable(subj_opt_vars,sprintf('all_fitted_subj_data_%s.csv',date));
        end
    end
end



