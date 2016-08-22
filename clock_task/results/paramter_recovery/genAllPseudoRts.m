uniform_flag = 1;
auto_corr_flag = 1;
% agent_names = {'fixedLR_softmax' 'kalman_uv_sum' 'kalman_softmax' 'qlearning'...
%     'franktc' 'kalman_processnoise' 'kalman_uv_logistic'};

agent_names = {'fixed_decay'};

for i = 1:length(agent_names)
    agent = agent_names{i};
    gen_temperatureFixed_rts(agent,uniform_flag,auto_corr_flag)
end