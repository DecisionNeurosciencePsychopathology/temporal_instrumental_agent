function s=clock_fit_model
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

[fittedparameters,options]=simps('clock_smooth_action_model',[.05 .95 100],[1 2 3],[options],[0.01 0.9 0.01],[0.3 .99 100],fargs{:});
alpha=fittedparameters(1); lambda=fittedparameters(2); epsilon=fittedparameters(3);

s.options = options;
s.parameters = fittedparameters;
[cost, constr,value_all]=clock_smooth_action_model(fittedparameters);

s.alpha = alpha;
s.lambda = lambda;
s.epsilon = epsilon;
s.value_all = value_all;
s.cost = cost;
% 
% if min([cost1 cost2 cost3])==cost1
%     %   s.prob.gamma=gamma1;
%     s.prob.alphawin=alphawin1;
%     s.prob.alphaloss=alphaloss1;
%     %   s.prob.WinAlpha=WinAlpha1;
%     %  s.prob.LossAlpha=LossAlpha1;
%     s.prob.c=c1;
%     s.prob.cost=cost1;
%     s.prob.e1=e11;
%     s.prob.e2=e21;
%     s.prob.e3=e31;
%     s.prob.prob1=prob11;
%     s.prob.prob2=prob21;
%     s.prob.prob3=prob31;
%     s.prob.delta=delta1;
%     s.prob.echosen=echosen1;
%     s.prob.etotal=etotal1;
% elseif min([cost1 cost2 cost3])==cost2
%     %   s.prob.gamma=gamma2;
%     s.prob.alphawin=alphawin2;
%     s.prob.alphaloss=alphaloss2;
%     %   s.prob.WinAlpha=WinAlpha2;
%     % s.prob.LossAlpha=LossAlpha2;
%     s.prob.c=c2;
%     s.prob.cost=cost2;
%     s.prob.e1=e12;
%     s.prob.e2=e22;
%     s.prob.e3=e32;
%     s.prob.prob1=prob12;
%     s.prob.prob2=prob22;
%     s.prob.prob3=prob32;
%     s.prob.delta=delta2;
%     s.prob.echosen=echosen2;
%     s.prob.etotal=etotal2;
% elseif min([cost1 cost2 cost3])==cost3
%     %   s.prob.gamma=gamma2;
%     s.prob.alphawin=alphawin3;
%     s.prob.alphaloss=alphaloss3;
%     %   s.prob.WinAlpha=WinAlpha2;
%     % s.prob.LossAlpha=LossAlpha2;
%     s.prob.c=c3;
%     s.prob.cost=cost3;
%     s.prob.e1=e13;
%     s.prob.e2=e23;
%     s.prob.e3=e33;
%     s.prob.prob1=prob13;
%     s.prob.prob2=prob23;
%     s.prob.prob3=prob33;
%     s.prob.delta=delta3;
%     s.prob.echosen=echosen3;
%     s.prob.etotal=etotal3;
% end

fprintf('alpha=%3f, lambda=%3f, epsilon=%3f,cost=%3f\n', s.alpha, s.lambda, s.epsilon, s.cost);
%figure(3); clf;
%plot(s.value_all); 
% save(sprintf('rl%d',id),'-struct', 's');

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % basic idea in curve fitting
% % is to have a function which returns
% % a goodness given some parameters
