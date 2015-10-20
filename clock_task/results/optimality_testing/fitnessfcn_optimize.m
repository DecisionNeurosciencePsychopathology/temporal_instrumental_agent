function [cost] = fitnessfcn_optimize(params, agent, runarray)
    %Numerical mapping system for running GA on each agentNOW and conditionNOW
    %%%%agentNOW%%%%
    %1 = kalman UV
    %2 = kalman skeptic
    %3 = kalman GRW
    %4 = Q-Leaning
    %5 = SARSA
    %%%%conditionNOW%%%%
    %1 = IEV
    %2 = DEV 
    %3 = QUADUP
    %4 = ALL
    
    %global agentNOW %which agent is being optimized
    %global conditionNOW %which lookup table is being used
    %global vperm %load in the permuted matrix
    
    %m=conditionNOW; %grab lookup
    %nruns = 15;
    %seed = 888;
    %ntrials = 200;
    %nbasis = 24;
    %ntimesteps = 500;
    
    %load in needed matrices
    % load('clock_options.mat')
    
    %Set up params for TD learning
%     if agentNOW >=4
%         clock_options.gamma = params(1);
%         clock_options.alpha = params(2);
%         clock_options.epsilon = params(3);
%         clock_options.lambda = params(4);
%         clock_options.gridcols = ntimesteps/10;
%         clock_options.episodeCount = ntrials;
%     end
   

    if strcmpi(agent.name, 'skeptic') || strcmpi(agent.name, 'skeptic_uvsum')
        cost = multirun_clock_sceptic(params, agent, runarray);
    elseif strcmpi(agent.name, 'kalman_grw')
        
    elseif strcmpi(agent.name, 'qlearning')
        
    elseif strcmpi(agent.name, 'sarsa')
        
    elseif strcmpi(agent.name, 'frank_tc')
        
    else
        error('I dont'' recognize this agent');
    end
        
    
    %Meat and potatoes i.e. optimize specific model for specific condition
%     if isstruct(m)
%         if agentNOW ==1
%             
%         elseif agentNOW==2
%             cost = multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m.name, ntrials, nbasis, ntimesteps, 0, 0, m},vperm);
%         elseif agentNOW==3
%             cost = multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m.name, ntrials, nbasis, ntimesteps, 0, 0, m},vperm);
%         elseif agentNOW==4
%             clock_options.agent = 'qlearning';
%             cost = multirun_cliff_walker_optimize(nruns, seed, {clock_options,m},vperm);
%         else
%             clock_options.agent = 'sarsa';
%             cost = multirun_cliff_walker_optimize(nruns, seed, {clock_options,m},vperm);
%         end
%     else %Whenever conditionNOW is set to 'all'
%         
%         %I hate breaking the DRY principle, and I probably could just make
%         %the above portion a modified for loop, but this works for now 
%         
%         cost=0; %Initialize cost
%         
%         %Run through each condition in the same run to get optimal params for 'all' conditions
%         if agentNOW ==1
%             
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 1, m{1}},vperm);
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{2}.name, ntrials, nbasis, ntimesteps, 0, 1, m{2}},vperm);
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {params, 0, m{3}.name, ntrials, nbasis, ntimesteps, 0, 1, m{3}},vperm);
%         elseif agentNOW==2
%             
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 0, m{1}},vperm);
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m{2}.name, ntrials, nbasis, ntimesteps, 0, 0, m{2}},vperm);
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0 0.2], 0, m{3}.name, ntrials, nbasis, ntimesteps, 0, 0, m{3}},vperm);
%         elseif agentNOW==3
%             
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m{1}.name, ntrials, nbasis, ntimesteps, 0, 0, m{1}},vperm);
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m{2}.name, ntrials, nbasis, ntimesteps, 0, 0, m{2}},vperm);
%             cost = cost + multirun_clock_logistic_operator_kalman_optimize(nruns, seed, {[params 0.9 0.2], 0, m{3}.name, ntrials, nbasis, ntimesteps, 0, 0, m{3}},vperm);
%         elseif agentNOW==4
%             
%             clock_options.agent = 'qlearning';
%             cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{1}},vperm);
%             cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{2}},vperm);
%             cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{3}},vperm);
%         else
%             
%             clock_options.agent = 'sarsa';
%             cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{1}},vperm);
%             cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{2}},vperm);
%             cost = cost + multirun_cliff_walker_optimize(nruns, seed, {clock_options,m{3}},vperm);
%         end
%     end
end
