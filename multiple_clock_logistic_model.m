function s=multiple_clock_logistic_model
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

fmincon_options = optimset(@fmincon);

opts = optimset('fmincon');
opts.LargeScale = 'off';
opts.Algorithm = 'active-set';
opts.Display = 'none';

%                 [params, SE, exitflag, xstart] = rmsearch(@(params) TC_minSE(params, subjdata{f}, model, emoSubset), 'fmincon', init_params, ...
%                     lower_limits, upper_limits, 'initialsample', num_start_pts, 'options', opts);
num_start_pts=50;
init_params_1 = [.01 .9877 -.06];
lower_bounds = [0.001 0.9 -1];
upper_bounds = [0.2 .999 0];

% init_params_1 = [.01 -.06];
% lower_bounds = [0.001 -0.3];
% upper_bounds = [0.2 0];

init_params_1 = [-.06];
lower_bounds = [-0.4];
upper_bounds = [0];

% init_params_1 = [.1];
% lower_bounds = [0];
% upper_bounds = [1];

%[fittedparameters_rmsearch,cost_rmsearch,exitflag_rmsearch,xstart_rmsearch]=rmsearch(@(params) clock_logistic_operator(params), 'fmincon', init_params_1, lower_bounds, upper_bounds, 'initialsample', num_start_pts, 'options', opts);

epsvalues = -0.6:.003:-.01;
%alphavalues = 0:.005:1;
costs=[];
for i = 1:length(epsvalues)
    %costs(i) = clock_logistic_operator(alphavalues(i));
    costs(i) = clock_logistic_operator(epsvalues(i));
end

%plot(alphavalues, costs);
%figure(2);
plot(epsvalues, costs);


[fittedparameters_fmincon, cost_fmincon, exitflag_fmincon] = fmincon(@(params) clock_logistic_operator(params), init_params_1, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

[fittedparameters_1,options]=simps('clock_logistic_operator', init_params_1, [1 2 3], options, lower_bounds, upper_bounds,fargs{:});
%[fittedparameters_1,options]=simps('clock_logistic_operator', init_params_1, [1 2], options, lower_bounds, upper_bounds,fargs{:});

alpha_1=fittedparameters_1(1); lambda_1=fittedparameters_1(2); epsilon_1=fittedparameters_1(3);
[cost_1, constr,value_all_1,value_hist_1]=clock_logistic_operator(fittedparameters_1);

fittedparameters_2=simps('clock_logistic_operator',[.1 .9 -.01],[1 2 3], options, [0.001 0.9 -1], [0.2 .99 0],fargs{:});
alpha_2=fittedparameters_2(1); lambda_2=fittedparameters_2(2); epsilon_2=fittedparameters_2(3);
[cost_2, constr,value_all_2,value_hist_2]=clock_logistic_operator(fittedparameters_2);

fittedparameters_3=simps('clock_logistic_operator',[.05 .99 -0.2],[1 2 3], options, [0.001 0.9 -1], [0.2 .99 0],fargs{:});
alpha_3=fittedparameters_3(1); lambda_3=fittedparameters_3(2); epsilon_3=fittedparameters_3(3);
[cost_3, constr,value_all_3,value_hist_3]=clock_logistic_operator(fittedparameters_3);

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
    s.value_all = value_all_1;
    s.cost = cost_1;
    %s.options = options;
    s.parameters = fittedparameters_1;
    s.value_hist=value_hist_1;
elseif min([cost_1 cost_2 cost_3])==cost_2
    s.alpha = alpha_2;
    s.lambda = lambda_2;
    s.epsilon = epsilon_2;
    s.value_all = value_all_2;
    s.cost = cost_2;
    %s.options = options;
    s.parameters = fittedparameters_2;
    s.value_hist=value_hist_2;
elseif min([cost_1 cost_2 cost_3])==cost_3
    s.alpha = alpha_3;
    s.lambda = lambda_3;
    s.epsilon = epsilon_3;
    s.value_all = value_all_3;
    s.cost = cost_3;
    %s.options = options;
    s.parameters = fittedparameters_3;
    s.value_hist=value_hist_3;
end

fprintf('alpha=%3f, lambda=%3f, epsilon=%3f,cost=%3f\n', s.alpha, s.lambda, s.epsilon, s.cost);
%figure(3); clf;
%plot(s.value_all);
% save(sprintf('rl%d',id),'-struct', 's');

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % basic idea in curve fitting
% % is to have a function which returns
% % a goodness given some parameters
