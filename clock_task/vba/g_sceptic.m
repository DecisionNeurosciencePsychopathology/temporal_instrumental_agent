function  [ gx ] = g_sceptic(x_t,phi,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT


beta = exp(phi(1));


gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;


v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

v_func = sum(v); %subjective value by timestep as a sum of all basis functions

p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature

rt_prev = u(1); %% retrieve previous RT

if strcmp(inG.autocorrelation,'exponential')
    lambda =  1./(1+exp(-phi(2))); %% introduce a choice autocorrelation parameter lambda
    chi =  1./(1+exp(-phi(3))); %% control the extent of choice autocorrelation

    p_choice = p_choice + chi.*(lambda.^(abs((1:ntimesteps) - rt_prev)));  %% incorporate an exponential choice autocorrelation function
    
    
    
%     rt_prev = 25;
%     lambdas = 0:.1:1;
%     for lambda = lambdas
%         lambda_val=lambda.^(abs((1:ntimesteps) - rt_prev));
%         plot(lambda_val)
%         hold on
%     end
    
    p_choice = p_choice./(sum(p_choice));  %% re-normalize choice probability so that it adds up to 1
elseif strcmp(inG.autocorrelation,'softmax_multitrial')
    
    p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature
end
    gx = p_choice';
end



