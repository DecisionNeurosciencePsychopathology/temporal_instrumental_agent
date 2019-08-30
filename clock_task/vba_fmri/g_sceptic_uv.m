function  [ gx ] = g_sceptic_uv(x_t, phi, u, inG)
% INPUT
% - x_t : hidden states (weights of basis functions)
% - phi : temperature (1x1)
% - u   : imput vector (not used in observation)
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t)

phi = transform_phi(phi, inG);
beta = phi(1); %force temperature to be positive through exponential
tau = phi(2); %for fixed_uv (not baked), tau operates in the observation function

gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;

hidden_state_index=1:inG.hidden_states*inG.nbasis; %total number of hidden states (inF.hidden_states is the number of state vectors)
hidden_state_index = reshape(hidden_state_index, inG.nbasis, inG.hidden_states); %3 x nbasis here

%V is first set of hidden states
v = x_t(hidden_state_index(:,1)) * ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

v_func = sum(v); %subjective value by timestep as a sum of all basis functions

%U is second set of hidden states
u = x_t(hidden_state_index(:,2)) * ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

u_func = sum(u); %subjective uncertainty by timestep as a sum of all basis functions

uv_func = v_func + tau.*u_func; %elementwise sum

p_choice = (exp((uv_func-max(uv_func))/beta)) / (sum(exp((uv_func-max(uv_func))/beta))); %Divide by temperature
gx = p_choice';

end
