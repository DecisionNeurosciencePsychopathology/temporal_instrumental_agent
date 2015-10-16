beta=.1;
model = 'v_processnoise';
ct = 1;
clear cost
for omega = 0.05:0.05:15
    [cost(ct),~]=skeptic_fitsubject_all_models_new_peSelect([.2261 beta omega],rts(1:50)*10,rews,seed, 24, 500, 0, 0, 500,model,1,2552);
    ct=ct+1;
end

figure(874); clf;
plot(0.05:0.05:omega,cost)




%See if giving the model the same sigma noise will return the omega used to
%generate the rts



%Plot cost as a function of sigma noise..
beta=.1;
model = 'v_processnoise';
ct = 1;
clear cost
omega=fitted_vars.noise_params.v_processnoise.best_parameters(1,:);
seed = [15 70 98 11];
rts = temperature_test_data.v_processnoise.temperature_rts.set_1;
rews = temperature_test_data.v_processnoise.temperature_rews.set_1;
for sigma_noise = 1:50:4000
    [cost(ct),~]=skeptic_fitsubject_all_models_new_peSelect([.2261 beta omega],rts(1:50)*10,rews,seed, 24, 500, 0, 0, 500,model,1,sigma_noise);
    ct=ct+1;
end
figure(854); clf;
plot(0:50:sigma_noise,cost)
hold on
plot(temperature_test_data.v_processnoise.ret.set_1.sigma_nosie(1),'rx')