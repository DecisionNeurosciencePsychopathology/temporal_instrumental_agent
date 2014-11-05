function s=multiple_wtw_logistic_model
load b;

% start off the learning
% algorithm with some initial values
% save modeldata.mat s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call Amoeba with our data
options(1)=1;        % To display intermediate results use 1, otherwise use 0 (default)
options(2)=1e-3;     % Relative x-tolerance
options(3)=1e-3;     % Relative f-tolerance
options(14)=100;    % Max. number of f-evaluations per internal
fargs={};

tic
% [fittedparameters_1,options]=simps('clock_logistic_operator',[.05 .95 1],[1 2 3],[options],[0.001 0.9 0.01],[0.2 .99 10],fargs{:});
% alpha_1=fittedparameters_1(1); lambda_1=fittedparameters_1(2); epsilon_1=fittedparameters_1(3);
% [cost_1, constr,value_all_1,value_hist_1]=clock_smooth_action_model(fittedparameters_1);

%epsilon is now scaled as the indifference point along mean(u_all) where agent switches from exploration to exploitation
%so should be a number between 0 and -1, roughly


%distrib_num = 1 - uniform, 2 - gen. Pareto, 3 - early beta, 4 - late beta
%5 - piecewise normal

distrib_num = 4;

[fittedparameters_1,options]=simps('wtw_logistic_operator',[0.01 0.9877 -2 -3],[1 2 3 4], options, [0.001 0.9 -6 -20], [0.2 .999 0 0],fargs{:},distrib_num);
alpha_1=fittedparameters_1(1); lambda_1=fittedparameters_1(2); epsilon_1=fittedparameters_1(3); k_1=fittedparameters_1(4);
[cost_1, constr,value_all_1,value_hist_1]=wtw_logistic_operator(fittedparameters_1,distrib_num);

fittedparameters_2=simps('wtw_logistic_operator',[.1 .9 -2 -7],[1 2 3 4], options, [0.001 0.9 -6 -20], [0.2 .99 0 0],fargs{:},distrib_num);
alpha_2=fittedparameters_2(1); lambda_2=fittedparameters_2(2); epsilon_2=fittedparameters_2(3); k_2=fittedparameters_2(4);
[cost_2, constr,value_all_2,value_hist_2]=wtw_logistic_operator(fittedparameters_2,distrib_num);

fittedparameters_3=simps('wtw_logistic_operator',[.05 .99 -2 -5],[1 2 3 4], options, [0.001 0.9 -6 -20], [0.2 .99 0 0],fargs{:},distrib_num);
alpha_3=fittedparameters_3(1); lambda_3=fittedparameters_3(2); epsilon_3=fittedparameters_3(3); k_3=fittedparameters_3(4);
[cost_3, constr,value_all_3,value_hist_3]=wtw_logistic_operator(fittedparameters_3,distrib_num);

toc
%
% s.options = options;
% s.parameters = fittedparameters;
% [cost, constr,value_all]=clock_smooth_action_model(fittedparameters);
%
% s.alpha = alpha;
% s.lambda = lambda;
% s.epsilon = epsilon;
% s.value_all = value_all;
% s.cost = cost;
%
if min([cost_1 cost_2 cost_3])==cost_1
    s.alpha = alpha_1;
    s.lambda = lambda_1;
    s.epsilon = epsilon_1;
    s.k = k_1;
    s.value_all = value_all_1;
    s.cost = cost_1;
    %s.options = options;
    s.parameters = fittedparameters_1;
    s.value_hist=value_hist_1;
elseif min([cost_1 cost_2 cost_3])==cost_2
    s.alpha = alpha_2;
    s.lambda = lambda_2;
    s.epsilon = epsilon_2;
    s.k = k_2;
    s.value_all = value_all_2;
    s.cost = cost_2;
    %s.options = options;
    s.parameters = fittedparameters_2;
    s.value_hist=value_hist_2;
elseif min([cost_1 cost_2 cost_3])==cost_3
    s.alpha = alpha_3;
    s.lambda = lambda_3;
    s.epsilon = epsilon_3;
    s.k = k_3;
    s.value_all = value_all_3;
    s.cost = cost_3;
    %s.options = options;
    s.parameters = fittedparameters_3;
    s.value_hist=value_hist_3;
end

fprintf('alpha=%3f, lambda=%3f, epsilon=%3f, k=%3f,cost=%3f\n', s.alpha, s.lambda, s.epsilon, s.k, s.cost);
%figure(3); clf;
%plot(s.value_all);
% save(sprintf('rl%d',id),'-struct', 's');

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % basic idea in curve fitting
% % is to have a function which returns
% % a goodness given some parameters
