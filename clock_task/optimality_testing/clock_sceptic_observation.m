function  [ gx ] = clock_sceptic_observation(x_t,P,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - u : [useless]
% - inG : 
% OUTPUT
% - gx : p(chosen|x_t)

beta = exp(P);
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;

v=x_t*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

v_func = sum(v); %subjective value by timestep as a sum of all basis functions
p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature

if inG.multinomial
gx = p_choice';
else
    gx = find(find(max(p_choice)),1);
end

