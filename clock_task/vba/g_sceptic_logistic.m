function  [ gx ] = g_sceptic_logistic(x_t,P,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - u : [useless]
% - inG :
% OUTPUT
% - gx : p(chosen|x_t) or RT

beta = exp(P(1));
%Need to add discrim as an additional param
discrim = 1./(1+exp(-P(2)));
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
u=x_t(nbasis+1:nbasis*2)*ones(1,ntimesteps) .* gaussmat; %Uncertainty is a function of Kalman uncertainties.

v_func = sum(v); %subjective value by timestep as a sum of all basis functions
u_func = sum(u); %vecotr of uncertainty by timestep


p_rt_exploit = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature
p_rt_explore = (exp((u_func-max(u_func))/beta)) / (sum(exp((u_func-max(u_func))/beta))); %Divide by temperature


%compared to other models that use a curve over which to choose (either by softmax or egreedy selection),
%kalman_uv_logistic computes explore and exploit choices and chooses according to a logistic.




u_final = sum(u_func)/length(u_func);

% if u_final == 0
%     rt_explore = ceil(.5*ntimesteps);
% else
%     rt_explore = find(u_func==max(u_func), 1); %return position of first max (and add gaussian noise?)
% end

sigmoid = 1/(1+exp(-discrim.*(u_final - inG.u_threshold))); %Rasch model with tradeoff as difficulty (location) parameter

p_choice_final = (((1 - sigmoid).*p_rt_explore) +  (sigmoid.*p_rt_exploit));

% if rand < sigmoid
%     %explore according to hardmax u
%     gx = rt_explore;
% else
%     gx = rt_exploit;
% end


if inG.multinomial
    gx = p_choice_final';
else
    best_rts = find(p_choice_final==max(p_choice_final));
    gx = mean(best_rts);
end