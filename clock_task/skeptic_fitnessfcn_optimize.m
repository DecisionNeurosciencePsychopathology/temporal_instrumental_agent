function [cost] = skeptic_fitnessfcn_optimize(params)

%hopefully this fixes the error of cell content reference
clear behav
global behavfiles sub model_names i

behav{sub}.data = readtable(behavfiles{sub},'Delimiter',',','ReadVariableNames',true);

%Initialize cost, seeds, and runs
cost=0;
seed = 66; %Let's aspire for greatness
nruns = 1; %Was 15 don't think thats really needed since the same cost is
%being computed

%NO
trialplots=0;

%Set model name for runs
modelname = model_names{i};

%9/9/15 JW: fix prop spread see what happens
params = [.2261 params];

%Sub loop code, does the heavy lifting
runs=unique(behav{sub}.data.run);
for run=1:length(runs);
cond = unique(behav{sub}.data.rewFunc(behav{sub}.data.run==runs(run)));
fprintf('\r%s subject#%d\r',char(cond),sub);
block = sprintf('%s_%d',char(cond),run);
rts_obs = behav{sub}.data.rt(behav{sub}.data.run==runs(run));
% apparently some RTs > 4000ms are recorded in the data: round them down to
% 4000ms
rts_obs(rts_obs>4000) = 4000;
rew_obs = behav{sub}.data.score(behav{sub}.data.run==runs(run));
% cond = behav{sub}.data.rewFunc(sub);
%% fit the model with a full range of RTs
% params = [.9165 .2261 .5]; %epsilon, prop_spread, spotlight


% if fit_params

    %clock_logistic_fitsubsubect(params, rts_obs', rew_obs');
    range = 'full';
    cost = cost + multirun_skeptic_fitsubsubect_all_models(nruns, seed, {params, rts_obs', rew_obs', 0, 24, 400, trialplots, 25, 400, modelname});
    sprintf('\n',params);
    %For ref delete later
    %cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 1, m{1}},vperm);      
    
    %% write p_choice fit statistics
    %behav{sub}.data.(sprintf('%s_p_chosen',modelname))(behav{sub}.data.run==runs(run),1) = ret.p_chosen';
end










% 
% %Numerical mapping system for running GA on each agentNOW and conditionNOW
% %%%%agentNOW%%%%
% %1 = kalman UV
% %2 = kalman skeptic
% %3 = kalman GRW
% %4 = Q-Leaning
% %5 = SARSA
% %%%%conditionNOW%%%%
% %1 = IEV
% %2 = DEV
% %3 = QUADUP
% %4 = ALL
% 
% global agentNOW %which agent is being optimized
% global conditionNOW %which lookup table is being used
% global vperm %load in the permuted matrix
% 
% m=conditionNOW; %grab lookup
% nruns = 15;
% seed = 888;
% ntrials = 200;
% nbasis = 24;
% ntimesteps = 500;
% 
% %load in needed matrices
% load('clock_options.mat')
% 
% %Set up params for TD learning
% if agentNOW >=4
%     clock_options.gamma = params(1);
%     clock_options.alpha = params(2);
%     clock_options.epsilon = params(3);
%     clock_options.lambda = params(4);
%     clock_options.gridcols = ntimesteps/10;
%     clock_options.episodeCount = ntrials;
% end
% 
% 
% %Meat and potatoes i.e. optimize specific model for specific condition
% if isstruct(m)
%     if agentNOW ==1
%         cost = multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m.name, ntrials, nbasis, ntimesteps, 0, 1, m},vperm);
%     elseif agentNOW==2
%         cost = multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m.name, ntrials, nbasis, ntimesteps, 0, 0, m},vperm);
%     elseif agentNOW==3
%         cost = multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m.name, ntrials, nbasis, ntimesteps, 0, 0, m},vperm);
%     elseif agentNOW==4
%         clock_options.agent = 'qlearning';
%         cost = multirun_cliff_walker_optimize(nruns, seed, {clock_options,m},vperm);
%     else
%         clock_options.agent = 'sarsa';
%         cost = multirun_cliff_walker_optimize(nruns, seed, {clock_options,m},vperm);
%     end
% else %Whenever conditionNOW is set to 'all'
%     
%     %I hate breaking the DRY principle, and I probably could subust make
%     %the above portion a modified for loop, but this works for now
%     
%     cost=0; %Initialize cost
%     
%     %Run through each condition in the same run to get optimal params for 'all' conditions
%     if agentNOW ==1
%         
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 1, m{1}},vperm);
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{2}.name, ntrials, nbasis, ntimesteps, 0, 1, m{2}},vperm);
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{3}.name, ntrials, nbasis, ntimesteps, 0, 1, m{3}},vperm);
%     elseif agentNOW==2
%         
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 0, m{1}},vperm);
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m{2}.name, ntrials, nbasis, ntimesteps, 0, 0, m{2}},vperm);
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m{3}.name, ntrials, nbasis, ntimesteps, 0, 0, m{3}},vperm);
%     elseif agentNOW==3
%         
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 0, m{1}},vperm);
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m{2}.name, ntrials, nbasis, ntimesteps, 0, 0, m{2}},vperm);
%         cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m{3}.name, ntrials, nbasis, ntimesteps, 0, 0, m{3}},vperm);
%     elseif agentNOW==4
%         
%         clock_options.agent = 'qlearning';
%         cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{1}},vperm);
%         cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{2}},vperm);
%         cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{3}},vperm);
%     else
%         
%         clock_options.agent = 'sarsa';
%         cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{1}},vperm);
%         cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{2}},vperm);
%         cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{3}},vperm);
%     end
% end
% end
