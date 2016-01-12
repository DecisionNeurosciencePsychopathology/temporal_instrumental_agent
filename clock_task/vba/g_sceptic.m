function  [ gx ] = g_sceptic(x_t,P,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - u : [useless]
% - inG : 
% OUTPUT
% - gx : p(chosen|x_t) or RT 

beta = exp(P);
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;

v=x_t*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

v_func = sum(v); %subjective value by timestep as a sum of all basis functions
p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature

if inG.multinomial
gx = p_choice';
else
    best_rts = find(p_choice==max(p_choice));
    gx = mean(best_rts);
end

