uniform_flag = 1;
auto_corr_flag = 0;
niv_flag=1;
% agent_names = {'fixedLR_softmax' 'kalman_uv_sum' 'kalman_softmax' 'qlearning'...
%     'franktc' 'kalman_processnoise' 'kalman_uv_logistic'};

%agent_names = {'fixed_decay'};


 agent_names = {'fixedLR_softmax' 'fixed_uv' 'fixed_decay' 'kalman_uv_sum' ...
       'kalman_softmax' 'kalman_processnoise' 'kalman_uv_logistic'};

for i = 1:length(agent_names)
    agent = agent_names{i};
    gen_temperatureFixed_rts(agent,uniform_flag,auto_corr_flag,niv_flag)
end