function  [ gx ] = g_sceptic(x_t,phi,u,inG)
% INPUT
% - x_t : 
% - phi : temperature (1x1)
% - u: not used
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t)

beta = exp(phi(1));
tau = phi(2);

gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

hidden_state_index=1:inG.hidden_state*nbasis; %total number of hidden states (inF.hidden_state is the number of state vectors)
hidden_state_index = reshape(hidden_state_index,nbasis,inG.hidden_state); %3 x nbasis here

mu=x_t(hidden_state_index(:,1)); %Value state vector
sigma = x_t(hidden_state_index(:,2)); %Uncertainty state vector

v=mu*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
v_func = sum(v); %subjective value by timestep as a sum of all basis functions

u=sigma*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
u_func = sum(u); %subjective uncertainty by timestep as a sum of all basis functions

uv_func = v_func + tau.*u_func;

p_choice = (exp((uv_func-max(uv_func))/beta)) / (sum(exp((uv_func-max(uv_func))/beta))); %Divide by temperature
gx = p_choice';
disp(gx)
%fprintf('gx: %s\n', num2str(gx));
if any(isnan(gx)), disp(gx); end
if any(isinf(gx)), disp(gx); end
end
