rng;
%Load in prev results so we don't overwrite
load('temperature_test_data');
data_sets = 100;
cost_og = zeros(1,data_sets);
model = 'v_processnoise';
beta=.1; %.1 is a good estimate
switch model
    case 'fixed'
        params = .2;
    case 'fixed_rho'
        params = [.2 .05];
    case 'fixed_KL'
        params = [.2 .002 .0007];
    case 'v_processnoise'
        params = 10;
    case 'uv'
        params = .45;
end

%10-15-15 So we've decided that each contingency will have a fixed sigma
%noise value, this means we will have to generate the same reward vecotr
%for a condition to have the same std.

condition = 'IEV';

switch condition
    case 'IEV'
        cond_seed = 15;
end


for i = 1:data_sets
seed=[cond_seed randi(200,1,3)];

[cost_og(i), ret.(['set_' num2str(i)])] = skeptic_fitsubject_all_models_new_peSelect_test_model([.2261 beta params],(ones(1,50))',100,seed, 24, 500, 0, 0, 500,model);

temperature_test_data.(model).temperature_rts.(['set_' num2str(i)]) = ret.(['set_' num2str(i)]).rt_obs(1:50);
temperature_test_data.(model).temperature_rews.(['set_' num2str(i)]) = ret.(['set_' num2str(i)]).rew_obs;
%seeds(i,:) = seed; %this causes the same seed to be produced everytime
%what gives?
end

%Save costs and ret to to boot
temperature_test_data.(model).ret = ret;
temperature_test_data.(model).cost = cost_og;
%temperature_test_data.(model).seeds = seeds;
 

save temperature_test_data temperature_test_data


%See what some generated rts are
%[cost,ret]=skeptic_fitsubject_all_models_new_peSelect_test_model([.2261 beta params],(ones(1,50))',(zeros(1,50))',seed, 24, 500, 1, 0, 500,model);

%Plot some gen rts if you wish
% plot(temperature_test_data.(model).temperature_rts.set_11);
% hold on
% 
% plot(temperature_test_data.(model).temperature_rts.set_42);
% plot(temperature_test_data.(model).temperature_rts.set_77);