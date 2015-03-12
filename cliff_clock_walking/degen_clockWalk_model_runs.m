close all

%params(1)=-.09; %epsilon
params(1)=-1; %epsilon
params(2)=0.2; %alpha
params(3)=0.99; %lambda
rlModel_cost_outputs.logistic_params=params;
seeds=[98 83];
condy = {'IEV','DEV','QUADUP'};
trials=[50 200];
for j = 1:length(condy)
    condition=condy{j};
    for i = 1:length(trials)
        
        %Delta
        clock_options_generator(condition)
        load('clock_options.mat')
        options.episodeCount = trials(i);
        options.gridrows = 1; %Agent can only wait or quit
        options.goal.row = 1;
        %options.epsilon=1;
        %options.decayflag=1;
        options.waitflag=1;
        [deggen_cost]=ClockWalking_noTempGen(options);
        
        rlModel_cost_outputs.(['deltaRule_' num2str(options.episodeCount) '_' condition])=max(deggen_cost);
        h(1) = figure(10000);
        save_fig(h(1),['deltaRule_' num2str(options.episodeCount) '_' condition],'fig')
        close all
        
%         %Q
%         clock_options_generator(condition)
%         load('clock_options.mat')
%         options.episodeCount = trials(i);
%         [clockWalker_cost]=ClockWalking(options);
%         rlModel_cost_outputs.(['qLearning_op_' num2str(options.episodeCount) '_' condition])=max(clockWalker_cost);
%         h(2) = figure(10002);
%         save_fig(h(2),['qLearning_' num2str(options.episodeCount) '_' condition],'fig')
%         close all
        
        
        %SARSA
        clock_options_generator(condition)
        load('clock_options.mat')
        options.episodeCount = trials(i);
        options.agent = 'sarsa';
        %options.epsilon=1;
        %options.decayflag=1;
        options.waitflag=1;
        [clockWalker_cost2]=ClockWalking(options);
        rlModel_cost_outputs.(['sarsa_op_' num2str(options.episodeCount) '_' condition])=max(clockWalker_cost2);
        h(3) = figure(10002);
        save_fig(h(3),['SARSA_' num2str(options.episodeCount) '_' condition],'fig')
        close all
        
%         %logistic
%         [cost_smooth_op1]=clock_logistic_operator(params,seeds,condition,trials(i));
%         rlModel_cost_outputs.(['logisitic_op_' num2str(trials(i)) '_' cond])=cost_smooth_op1;
%         h(4) = figure(1);
%         save_fig(h(4),['logisticOP_' num2str(options.episodeCount) '_' condition],'fig')
%         close all
        
    end
end
    save rlModel_cost_outputs rlModel_cost_outputs