function gen_temperatureFixed_rts(agent)
%This script will generate a set of rts and rews to be used in future
%parameter recovery simulations.
%
% ex: gen_temperatureFixed_rts('v_processnoise')
%

%Load in prev results so we don't overwrite
load('C:\kod\temporal_instrumental_agent\clock_task\results\paramter_recovery\param_recovery_test_data');
%Load in master sample array
load('C:\kod\temporal_instrumental_agent\clock_task\results\optimality_testing\mastersamp.mat')
%Load in clock options
load('C:\kod\temporal_instrumental_agent\clock_task\optimality_testing\clock_options.mat')


data_sets = 100; %Define number of data sets
cost_og = zeros(1,data_sets);
%agent = 'kalman_sigmavolatility'; %Which agent to use
beta=.1; %.1 is a good estimate

%Set some parameters
rngseeds=[98 83 66 10];
ntrials=50;
nbasis=24;
ntimesteps=500;
reversal=0;

%Set agent parameters
switch agent
    case 'fixedLR_softmax'
        params = [0.03081894 1.536380918 0.01698200]; %prop_spread  beta  alpha
%     case 'fixed'
%         params = .2;
%     case 'fixed_rho'
%         params = [.2 .05];
%     case 'fixed_KL'
%         params = [.2 .002 .0007];
%     case 'v_processnoise'
%         params = 10;
%     case 'uv'
%         params = .45;
    case 'kalman_sigmavolatility'
        params = [0.272621000 0.10000000 0.00000000 0.3532713]; %All condition params via Michael's email
        %params = [0.272621000 0.10000000 0 0]; %All condition params via recovery
    case 'frank'
        params = [0.03059226 87511.530 3.449525 2.092848  685.891054 1879.997  611.3465]; %All condition params via Michael's email [lambda   epsilon   alphaG   alphaN  K   nu  rho]
        priors.V=0; %initialize expected value for first trial to prior (possibly from previous run)
        priors.Go=0; %initialize Go for first trial
        priors.NoGo=0; %initialize NoGo for first trial
        rtbounds = [1 5000]; %Don't let the agent choose 0
    case 'kalman_uv_sum'
        %All condition params via Michael's email
        params =[0.62416337 0.001524063 0.69507251]; %prop_spread  beta  discrim
    case 'kalman_uv_sum_kl'
        %All condition params via Michael's email
        params =[0.27298721 1.2574299 0.611199892 2.15673530 2.2536939758]; %prop_spread  beta  tau kappa lambda
    case 'kalman_uv_logistic'
        %All condition params via Michael's email
        params =[0.62416337 0.001524063 0.69507251]; %prop_spread  beta  tau
    case 'qlearning'
        %All condition params via Michael's email
        params =[]; %gamma, alpha, epsilon lambda
    otherwise
        error('Not any agent I''ve heard of');
end

%10-15-15 So we've decided that each contingency will have a fixed sigma
%noise value, this means we will have to generate the same reward vecotr
%for a condition to have the same std.

condition = 'IEV';
cond = mastersamp.(condition);

%Set repeatable range
rng(15);
condnames=fieldnames(mastersamp); %by default, optimize over all contingencies in lookup
ncond=length(condnames); %how many conditions/contingencies

%truncate mastersamp to the number of trials used here to save RAM
for i = 1:length(condnames)
    mastersamp.(condnames{i}).lookup = mastersamp.(condnames{i}).lookup(:,1:ntrials);
end

%Permuste the master lookup
    clear row
    for j = 1:data_sets
        tmp = mastersamp.(condition);
        tmp.name = condition; %for keeping track of what contingency this is (esp. when they mix in ALL)
        tmp.lookup = tmp.lookup(:, randperm(size(tmp.lookup,2))); %randomly permute columns
        %tmp.perm = randperm(size(tmp.lookup,2)); %just for checking that this works
        %tmp.lookup = tmp.lookup(:, tmp.perm); %randomly permute columns
        row(j) = tmp; %add to a vector of structs
    end
    optmat = row; 

    
    %Using mike's new optimality script
    for i = 1:data_sets
        if (strcmp(agent,'frank'))
            [cost_og(i),~, ret.(['set_' num2str(i)])]=TC_Alg_forward(params, priors, optmat(i), rngseeds, ntrials,rtbounds); %No rtbounds arg default is 0-5000ms
            ret.(['set_' num2str(i)]).rts = ret.(['set_' num2str(i)]).rtpred';
            ret.(['set_' num2str(i)]).rew_i = ret.(['set_' num2str(i)]).rew;
            ret.(['set_' num2str(i)]).optmat = optmat(i);
            [cost_og(i),~,~,~,ret.(['set_' num2str(i)])] = clock_sceptic_agent_genRts(params, agent, rngseeds, optmat(i), ntrials, nbasis, ntimesteps, reversal);
        elseif(strcmp(agent,'qlearning') || (strcmp(agent,'sarsa'))) 
            [cost_og(i),~,~,~,~, ret.(['set_' num2str(i)])]=ClockWalking_3D_discountedEv_genRts(clock_options,optmat,rngseeds,params);
        else
            [cost_og(i),~,~,~,ret.(['set_' num2str(i)])] = clock_sceptic_agent_genRts(params, agent, rngseeds, optmat(i), ntrials, nbasis, ntimesteps, reversal);
        end
    end


%Let use a new data structure to save results
param_recovery_test_data.(agent).ret = ret;
param_recovery_test_data.(agent).cost = cost_og;
param_recovery_test_data.(agent).condition = condition;
%param_recovery_test_data.(agent).optmat = optmat;


save('C:\kod\temporal_instrumental_agent\clock_task\results\paramter_recovery\param_recovery_test_data', 'param_recovery_test_data')


%See what some generated rts are
%[cost,ret]=skeptic_fitsubject_all_models_new_peSelect_test_model([.2261 beta params],(ones(1,50))',(zeros(1,50))',seed, 24, 500, 1, 0, 500,model);

%Plot some gen rts if you wish
% plot(temperature_test_data.(model).temperature_rts.set_11);
% hold on
% 
% plot(temperature_test_data.(model).temperature_rts.set_42);
% plot(temperature_test_data.(model).temperature_rts.set_77);


%Old code
% for i = 1:data_sets
% seed=[cond_seed randi(200,1,3)];
% 
% [cost_og(i), ret.(['set_' num2str(i)])] = skeptic_fitsubject_all_models_new_peSelect_test_model([.2261 beta params],(ones(1,50))',100,seed, 24, 500, 0, 0, 500,model);
% 
% temperature_test_data.(model).temperature_rts.(['set_' num2str(i)]) = ret.(['set_' num2str(i)]).rt_obs(1:50);
% temperature_test_data.(model).temperature_rews.(['set_' num2str(i)]) = ret.(['set_' num2str(i)]).rew_obs;
% %seeds(i,:) = seed; %this causes the same seed to be produced everytime
% %what gives?
% end

%Save costs and ret to to boot
% temperature_test_data.(model).ret = ret;
% temperature_test_data.(model).cost = cost_og;
%temperature_test_data.(model).seeds = seeds;