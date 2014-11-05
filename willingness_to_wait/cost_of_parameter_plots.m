%% exploration
eps=-10:.1:0;

for i = 1:length(eps)
    disp(eps(i))
    cost(i) = wtw_logistic_operator([0.18073 0.98811 eps(i) -13.957],4);
end
figure(2);clf;
plot(eps, cost);


%% learning rate
alpha=0:.005:.5;

for i = 1:length(alpha)
    disp(alpha(i))
    cost(i) = wtw_logistic_operator([alpha(i) 0.96058 -2.1033 -15.967],3);
end
figure(2);clf;
plot(alpha, cost);

%% discount rate
clear;
log_k=-1000:1:0;

for i = 1:length(log_k)
    disp(log_k(i))
    cost(i) = wtw_logistic_operator([0.18073 0.98811 -5.9234 log_k(i)],4);
end
figure(2);clf;
plot(log_k, cost);

%% lambda
clear;
lambda=.98:.0001:1;

for i = 1:length(lambda)
    disp(lambda(i))
    cost(i) = wtw_logistic_operator([0.18073 lambda(i) -5.9234 -13.957],4);
end
figure(2);clf;
plot(lambda, cost);