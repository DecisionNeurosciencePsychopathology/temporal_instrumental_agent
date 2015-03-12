function s=ClockWalker_fit_model
%load clock_options %??

% start off the learning
% algorithm with some initial values
% save modeldata.mat s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call Amoeba with our data
fit_options(1)=1;        % To display intermediate results use 1, otherwise use 0 (default)
fit_options(2)=1e-3;     % Relative x-tolerance
fit_options(3)=1e-3;     % Relative f-tolerance
fit_options(14)=100;    % Max. number of f-evaluations per internal
fargs={};

[fittedparameters,fit_options]=simps('ClockWalking_3D',[.05 .9 .2],[1 2 3],[fit_options],[0 0.9 0],[0.5 1 1],fargs{:});
alpha1=fittedparameters(1); gamma1=fittedparameters(2); epsilon1=fittedparameters(3);
[cost1, constr,quits1,cumReward1]=ClockWalking_3D(fittedparameters);

fittedparameters=simps('ClockWalking_3D',[.1 .99 .3],[1 2 3],[fit_options],[0 0.9 0],[0.5 1 1],fargs{:});
alpha2=fittedparameters(1); gamma2=fittedparameters(2); epsilon2=fittedparameters(3);
[cost2, constr,quits2,cumReward2]=ClockWalking_3D(fittedparameters);

[fittedparameters]=simps('ClockWalking_3D',[.2 .999 .4],[1 2 3],[fit_options],[0 0.9 0],[0.5 1 1],fargs{:});
alpha3=fittedparameters(1); gamma3=fittedparameters(2); epsilon3=fittedparameters(3);
[cost3, constr,quits3,cumReward3]=ClockWalking_3D(fittedparameters);

[fittedparameters]=simps('ClockWalking_3D',[.3 .9999 .5],[1 2 3],[fit_options],[0 0.9 0],[0.5 1 1],fargs{:});
alpha4=fittedparameters(1); gamma4=fittedparameters(2); epsilon4=fittedparameters(3);
[cost4, constr,quits4,cumReward4]=ClockWalking_3D(fittedparameters);


s.options = fit_options;
% s.parameters = fittedparameters;
% 
% 
% s.alpha = alpha;
% s.gamma = gamma;
% s.epsilon = epsilon;
% s.value_all = value_all;
% s.cost = cost;

alphas = [alpha1 alpha2 alpha3 alpha4];
gammas = [gamma1 gamma2 gamma3 gamma4];
epsilons = [epsilon1 epsilon2 epsilon3 epsilon4];


costs = [cost1 cost2 cost3 cost4];
min_idx = find(costs==min(costs));

s.cost = costs(min_idx);
s.alpha = alphas(min_idx);
s.gamma = gammas(min_idx);
s.epsilon = epsilons(min_idx);


fprintf('alpha=%3f, lambda=%3f, epsilon=%3f,cost=%3f\n', s.alpha, s.gamma, s.epsilon, s.cost);
%figure(3); clf;
%plot(s.value_all); 
% save(sprintf('rl%d',id),'-struct', 's');

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % basic idea in curve fitting
% % is to have a function which returns
% % a goodness given some parameters
