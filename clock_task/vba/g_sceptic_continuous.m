function  [ gx ] = g_sceptic_continuous(x_t,phi,u,inG)
% INPUT
% - x_t : value (nbasisx1)
% - phi(1) : K (1x1)
% - phi(2) : lambda (1x1)
% - u(3)   : previous RT (n_trialsx1
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : RT

K = unifinv(fastnormcdf(phi(1)), - inG.maxRT./2, inG.maxRT./2); %uniform transform
lambda = 1 ./ (1+exp(-phi(2))); %exponential transform to 0..1

gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

RT_prev = u(3);

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

v_func = sum(v); %subjective value by timestep as a sum of all basis functions

best_rts = mean(find(v_func==max(v_func)));

gx = wmean([best_rts + K, RT_prev], [1-lambda, lambda],2);


