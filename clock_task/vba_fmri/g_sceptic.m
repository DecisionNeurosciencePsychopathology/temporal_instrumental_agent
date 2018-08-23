function  [ gx ] = g_sceptic(x_t, phi, u, inG)
% INPUT
% - x_t : hidden states (weights of basis functions)
% - phi : temperature (1x1)
% - u   : imput vector (not used in observation)
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t)

phi = transform_phi(phi, inG);
beta = phi(1); %force temperature to be positive through exponential

gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

v_func = sum(v); %subjective value by timestep as a sum of all basis functions

p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature
gx = p_choice';

end
