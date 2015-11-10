function [cost] = fitnessfcn_optimize(params, agent, runarray)
    if ismember(agent.name, {'fixedLR_softmax', 'fixedLR_egreedy', 'fixedLR_egreedy_grw', ...
            'asymfixedLR_softmax', 'kalman_softmax', 'kalman_processnoise', ...
            'kalman_sigmavolatility', 'kalman_uv_logistic', 'kalman_uv_sum', ...
            'fixedLR_kl_softmax', 'kalman_kl_softmax', 'kalman_processnoise_kl', 'kalman_uv_sum_kl'})
        cost = multirun_clock_sceptic(params, agent, runarray);        
    elseif strcmpi(agent.name, 'qlearning')
        agent.clock_options.episodeCount = agent.ntrials;
        agent.clock_options.ntimesteps = agent.timesteps/10; %TD models operate on 50 timesteps, not 500, typically
        cost = multirun_cliff_walker_optimize(params, agent, runarray);
    elseif strcmpi(agent.name, 'sarsa')
        agent.clock_options.episodeCount = agent.ntrials;
        agent.clock_options.ntimesteps = agent.timesteps/10;
        cost = multirun_cliff_walker_optimize(params, agent, runarray);
    elseif strcmpi(agent.name, 'franktc')
        cost = multirun_TC_alg_forward(params, agent, runarray);
    else
        error('I dont'' recognize this agent');
    end
end
