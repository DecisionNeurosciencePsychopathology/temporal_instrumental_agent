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

%fmincon_options = optimset(@fmincon);
%fmincon_options = optimoptions(@fmincon, 'UseParallel',true, 'Algorithm', 'active-set', 'DiffMinChange', 0.01);
fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set', 'DiffMinChange', 0.005);

opts = optimset('fmincon');
opts.LargeScale = 'off';
opts.Algorithm = 'active-set';
opts.Display = 'none';

num_start_pts=50;
init_params_1 = [.01 .9877 -.06];
lower_bounds = [0.001 0.9 -1];
upper_bounds = [0.2 .999 0];

%alpha and epsilon only
% init_params_1 = [.01 -.06];
% lower_bounds = [0.001 -0.3];
% upper_bounds = [0.2 0];

%epsilon only
init_params_1 = [-.08];
init_params_2 = [-.08];
lower_bounds = [-0.5];
upper_bounds = [-.005];

% init_params_1 = [.1];
% lower_bounds = [0];
% upper_bounds = [1];

%[fittedparameters_rmsearch,cost_rmsearch,exitflag_rmsearch,xstart_rmsearch]=rmsearch(@(params) clock_logistic_operator(params), 'fmincon', init_params_1, lower_bounds, upper_bounds, 'initialsample', num_start_pts, 'options', opts);

poolobj = parpool('local', 4);
epsvalues = -0.6:.005:-.001;
%alphavalues = 0:.005:1;
costs=zeros(1,length(epsvalues));
rts=zeros(length(epsvalues), 125); %125 trials
for i = 1:length(epsvalues)
    %costs(i) = clock_logistic_operator(alphavalues(i));
    %[costs(i) dummy1 dummy2 dummy3 rts(i,:)] = clock_logistic_operator(epsvalues(i));
    [costs(i)] = multirun_clock_logistic_operator(20, 999, {epsvalues(i), 0, 'IEV', 100}); %clock_logistic_operator(epsvalues(i));
end
delete(poolobj);

%plot(alphavalues, costs);
figure(5);
plot(epsvalues, costs);

diff(costs)

figure(1);
plot(epsvalues, costs_shuffle);

sum(diff(costs_shuffle)==0)

figure(2);
plot(epsvalues, costs_hard);

sum(diff(costs_hard)==0)

figure(3);
plot(epsvalues, costs_soft);

sum(diff(costs_soft)==0)




%identical costs at epsilon = -.45 and -.455
[cost_1, dummy1, dummy2, dummy3, rts_1] = clock_logistic_operator(-.45);
[cost_2, dummy1, dummy2, dummy3, rts_2] = clock_logistic_operator(-.4550);


[cost_1, dummy1, dummy2, dummy3, rts_1] = clock_logistic_operator(-.08);


%testing rbfeval function
weights=[0.0776 7.4801 0.0792 0.0008 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0024 0.2286 22.0326 0.2332];
centers=[-454.5455 0.0000 454.5455 909.0909 1363.6364 1818.1818 2272.7273 2727.2727 3181.8182 3636.3636 4090.9091 4545.4545 5000.0000 5454.5455];
widths=ones(1,14).*350.0700;
 
 outs=[];
 for v = 0:5000
     outs(v+1) = rbfeval(v, weights, centers, widths);
end

plot(0:5000, outs);

% 
% calc=fmincon(@(params) rbfeval(params, weights, centers, widths), 1000, [], [], [], [], 0, 5000, [], opts);
% calc=fminbnd(@(params) -rbfeval(params, weights, centers, widths), 0, 5000);
% rbfeval(calc, weights, centers, widths)

%optionsGA = gaoptimset('PlotFcns',@gaplotbestfun,'PlotInterval',5, ...
%    'PopInitRange',[-5;5]);

optionsGA = gaoptimset('PlotFcns',@gaplotbestf,'Display','iter','InitialPopulation',init_params_1);
% We run GA with the options 'optionsGA' as the tenth argument.
numVars = 3;

init_params_1 = [-.06 .01 .9877 ];
lower_bounds = [-1 0.001 0.9 ];
upper_bounds = [0 0.2 .999];

[Xga,Fga] = ga(@fitnessfcn, ...
    numVars, ... %number of free parameters
    [], ... A*x
    [], ... %b
    [], ... %Aeq
    [], ... %beq
    lower_bounds, ... %LB
    upper_bounds, ... %UB
    [], ... %nonlcon
    optionsGA);

%best values from genetic algorithm:  -0.0600    0.0100    0.9875
gabest=[-0.0600    0.0100    0.9875];


%more best values from 20 runs of data: seed=999
gabest = [-0.0434    0.0202    0.9706];

[cost,constr,value_all,value_hist,rts,mov] = clock_logistic_operator(gabest, [20 80], 'IEV', 100);


%encourage optimizer to look around!
%fmincon_options.
%wrapper for multiple draws of clock_logistic_operator
poolobj = parpool('local', 4);
[test, totcost, costs] = fmincon(@(params) multirun_clock_logistic_operator(4, 999, params), init_params_1, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
delete(poolobj);

[totcost, runcosts, runseeds] = multirun_clock_logistic_operator(10, 999, init_params_1);

[fittedparameters_fmincon, cost_fmincon, exitflag_fmincon] = fmincon(@(params) clock_logistic_operator(params, [93 62], 'DEV', 125), init_params_1, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

%[fittedparameters_1,options]=simps('clock_logistic_operator', init_params_1, [1 2 3], options, lower_bounds, upper_bounds,fargs{:});
[fittedparameters_1,options]=simps('clock_logistic_operator', init_params_1, [1], options, lower_bounds, upper_bounds,fargs{:});

%alpha_1=fittedparameters_1(1); lambda_1=fittedparameters_1(2); epsilon_1=fittedparameters_1(3);
%[cost_1, constr,value_all_1,value_hist_1]=clock_logistic_operator(fittedparameters_1);

[fittedparameters_2,options]=simps('clock_logistic_operator', init_params_2, [1], options, lower_bounds, upper_bounds,fargs{:});
%fittedparameters_2=simps('clock_logistic_operator',[.1 .9 -.01],[1 2 3], options, [0.001 0.9 -1], [0.2 .99 0],fargs{:});
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


%% FRANK TC MODEL SECTION


    init_params = [ 0.3 ; 2000 ; 0.2 ; 0.2 ; 1000 ; 0.1 ; 0.5 ; 300 ];
    lower_limits = [ 0 ; 0 ; 0.01 ; 0.01 ; .1 ; 0 ; .1 ; 0 ];
    upper_limits = [1 ; 100000 ; 5 ; 5 ; 5000 ; 5000 ; 5000 ; 10000 ]; % for rmsearch set min/max to 0 for unused params (otherwise spits out weird values that aren't used)
    hdr = {'Subject','Session','lambda','explore','alphaG','alphaN','K','nu','ignore','rho','SSE'};

%lambda, epsilon, alphaG, alphaN, K, nu (go for gold), rho
%init_params = [ 0.3 ; 2000 ; 0.2 ; 0.2 ; 1000 ; 0.1 ; 300 ];
init_params = [ 0.2 ; 3000 ; 0.3 ; 0.3 ; 1000 ; 0.1 ; 300 ];
lower_bounds = [ 0 ; 0 ; 0.01 ; 0.01 ; 1 ; 0 ; 0 ];
upper_bounds = [1 ; 100000 ; 5 ; 5 ; 5000 ; 5000 ; 10000 ];

priors = [];
priors.V = 0;
priors.Go = 0;
priors.NoGo = 0;
[cost, RTpred, ret]=TC_Alg_forward(init_params, priors, 'DEV', 32, 500);


fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set', 'MaxFunEvals', 20000, 'TolFun', 1e-4, 'MaxIter', 2000);
[par, cost, exitflag] = fmincon(@(params) TC_Alg_forward(params, priors, 'DEV', 32, 100), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

fprintf('%.3f, ', par); fprintf('\n');


[recost, RTpred, ret]=TC_Alg_forward(par, priors, 'DEV', 32, 100);

figure(1); plot(RTpred);

%try same parameters for IEV
[recost, RTpred, ret]=TC_Alg_forward(par, priors, 'IEV', 32, 100);

figure(2); plot(RTpred);

(params, priors, cond, rngseeds, ntrials, rtbounds)

%see whether we can get reasonable params across multiple runs of each condition
%poolobj = parpool('local', 4);
[test, totcost, costs] = fmincon(@(params) multirun_TC_alg_forward({'IEV', 'DEV', 'IEV', 'QUADUP'}, 999, {params, priors, 'COND', 50}), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
%delete(poolobj);

fprintf('%.3f, ', test); fprintf('\n');


%from TC_Alg.m
% GENERATIVE model just pick some params to generate data
% global generative;
% if (generative ==1)
%     lambda = 0.2;
%     explore = 3000;
%     alpha1 = .3;
%     alpha2 = .3;
%     K = 1500;
%     scale = .25;
%     meandiff = 1000;
%     if strcmp(dist_type,'Gauss')
%         meandiff=20;
%         explore = 10;
%     end
%     Noise=2000;
% end

