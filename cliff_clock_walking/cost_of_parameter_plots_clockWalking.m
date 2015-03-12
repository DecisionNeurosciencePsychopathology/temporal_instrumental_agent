function cost_of_parameter_plots_clockWalking
%%Optimal values
load('clock_options.mat')
load('s.mat')
clock_op = initialize(clock_options,s); %Initalize optimal parms
clock_op.cond='QUADUP'; %set condition
clock_op.diag =0; %turn off plots


%% exploration
eps=0:.01:1;
cost = zeros(length(eps),1);
for i = 1:length(eps)
    disp(eps(i))
    clock_op.epsilon = eps(i);
    cost(i) = ClockWalking_3D(clock_op);
end
figure(2);clf;
plot(eps, cost);
title('Epsilon vs cost')
h(1)=figure(2);

%% learning rate
alpha=0:.01:.5;
cost = zeros(length(alpha),1);
clock_op = initialize(clock_op,s); %Initalize optimal parms
for i = 1:length(alpha)
    disp(alpha(i))
    clock_op.alpha = alpha(i);
    cost(i) = ClockWalking_3D(clock_op);
end
figure(3);clf;
plot(alpha, cost);
title('Alpha vs cost')
h(2)=figure(3);


% %% discount rate
% clear;
% log_k=-1000:1:0;
% 
% for i = 1:length(log_k)
%     disp(log_k(i))
%     cost(i) = ClockWalking_3D([0.18073 0.98811 -5.9234 log_k(i)],4);
% end
% figure(4);clf;
% plot(log_k, cost);

%% lambda

gamma=.9:.001:.999;
cost = zeros(length(gamma),1);
clock_op = initialize(clock_op,s); %Initalize optimal parms
for i = 1:length(gamma)
    disp(gamma(i))
    clock_op.gamma = gamma(i);
    cost(i) = ClockWalking_3D(clock_op);
end
figure(4);clf;
plot(gamma, cost);
title('Gamma vs cost')
h(3)=figure(4);


%Neat little save all figures trick
var_names = {'epsi', 'alpha', 'gamma'};
for i = 1:length(h)
save_fig(h(i),[clock_op.cond '_' var_names{i} '_cost_graph'],'fig')
end


function clock_op=initialize(options,opt_values)
options.alpha = opt_values.alpha;
options.gamma = opt_values.gamma;
options.epsilon = opt_values.epsilon;
clock_op=options;



