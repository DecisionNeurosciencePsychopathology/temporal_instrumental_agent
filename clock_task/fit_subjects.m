%script to load in subjects' data and fit using the logistic operator

cd('/Users/michael/Google_Drive/Data_Analysis/clock_analysis/fmri/behavior_files');
behavfiles = dir(fullfile(pwd, '*.csv'));

%read header
% fid = fopen(behavfiles(1).name, 'r');
% head = textscan(fid, '%s', 1);
% fclose(fid);
% m = csvread(behavfiles(1).name, 1, 0);

%wow, matlab has a useful function for mixed data types!!
data = readtable(behavfiles(1).name,'Delimiter',',','ReadVariableNames',true);

%just example first run
rts_obs = data.rt(1:50);
rew_obs = data.score(1:50);


params = [.02 .985 -.06]; %alpha, lambda, epsilon
clock_logistic_fitsubject(params, rts_obs', rew_obs');

fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set');%, 'DiffMinChange', 0.001);

init_params_1 = [.02 .9877 -.06];
lower_bounds = [0.001 0.9 -1];
upper_bounds = [0.2 .999 0];

[fittedparameters_fmincon, cost_fmincon, exitflag_fmincon] = fmincon(@(params) clock_logistic_fitsubject(params, rts_obs', rew_obs'), init_params_1, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

[~, fitted_object] = clock_logistic_fitsubject(-.09, rts_obs', rew_obs');
fitted_object.cost_total
[~, fitted_object] = clock_logistic_fitsubject(-.06, rts_obs', rew_obs');
fitted_object.cost_total

plot(1:50, fitted_object.rts_obs)
hold on;
plot(1:50, fitted_object.rts_pred, 'b')