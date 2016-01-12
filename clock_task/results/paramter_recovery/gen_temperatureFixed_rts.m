function gen_temperatureFixed_rts(agent,uniform_flag)
%This script will generate a set of rts and rews to be used in future
%parameter recovery simulations.
% IN:
%   agent -> name of agent to generate data
%   fname -> data sotrage matrix default is default_subj_test_data
%   ex: gen_temperatureFixed_rts('v_processnoise')
%

if nargin <2
    uniform_flag=0;
end

%Create the defualt data hub if it doesn't exisit
if ~exist('temporal_instrumental_agent\clock_task\results\paramter_recovery\param_recovery_test_data.mat','file')
    param_recovery_test_data = [];
    save('c:\kod\temporal_instrumental_agent\clock_task\results\paramter_recovery\param_recovery_test_data.mat','param_recovery_test_data'); %I think you need the full path so this may be different on your machine...
end




%Load in prev results so we don't overwrite
% try
%     load(['temporal_instrumental_agent\clock_task\results\paramter_recovery\',fname]);
% catch
%     fprintf('No pseudo subject data matrix found using default...\n\n');
%     load('temporal_instrumental_agent\clock_task\results\paramter_recovery\default_subj_test_data'); %This also maybe different depending on where default is saved
% end

%Load in prev results so we don't overwrite
load('temporal_instrumental_agent\clock_task\results\paramter_recovery\param_recovery_test_data');

%You will need these for certain models
%Load in master sample array
load('temporal_instrumental_agent\clock_task\optimality_testing\mastersamp.mat')
%Load in clock options
load('temporal_instrumental_agent\clock_task\optimality_testing\clock_options.mat')


data_sets = 100; %Define number of data sets
cost_og = zeros(1,data_sets);
%agent = 'kalman_sigmavolatility'; %Which agent to use
beta=.1; %.1 is a good estimate

%Set some parameters
rngseeds=[98 83 66 10];
ntrials=50;
clock_options.episodeCount=ntrials;
nbasis=24;
ntimesteps=500;
reversal=0;

a = initialize_stability_struct();

if uniform_flag
     [params,priors,rtbounds] = getAgentParamsUniform(agent,beta,data_sets,a);
else
    [params,priors,rtbounds] = getAgentParams(agent);
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
for i = 1:length(condnames)-1 %Quick dirty fix since DEVLINPROB qasn't correctly made
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
        
        %If there are multiple parameter values for each model loop through
        %them to create the rts
        if size(params,1)==1
            j=1;
        else
            j=i;
        end
        
        
        if (strcmp(agent,'franktc'))
            [cost_og(i),~, ret.(['set_' num2str(i)])]=TC_Alg_forward(params(j,:), priors, optmat(i), rngseeds, ntrials,rtbounds); %No rtbounds arg default is 0-5000ms
            ret.(['set_' num2str(i)]).rts = ret.(['set_' num2str(i)]).rtpred';
            ret.(['set_' num2str(i)]).rew_i = ret.(['set_' num2str(i)]).rew;
            ret.(['set_' num2str(i)]).optmat = optmat(i);
        elseif(strcmp(agent,'qlearning') || (strcmp(agent,'sarsa'))) 
            clock_options.agent=agent;
            [cost_og(i),~,~,~,~, ret.(['set_' num2str(i)])]=ClockWalking_3D_discountedEv_genRts(clock_options,optmat(i),rngseeds,params(j,:));
        else
            [cost_og(i),~,~,~,ret.(['set_' num2str(i)])] = clock_sceptic_agent_genRts(params(j,:), agent, rngseeds, optmat(i), ntrials, nbasis, ntimesteps, reversal);
        end
    end


%Let use a new data structure to save results
param_recovery_test_data.(agent).ret = ret;
param_recovery_test_data.(agent).cost = cost_og;
param_recovery_test_data.(agent).condition = condition;
param_recovery_test_data.(agent).params = params;
param_recovery_test_data.(agent).name = agent;
%param_recovery_test_data.(agent).optmat = optmat;


save('C:\kod\temporal_instrumental_agent\clock_task\results\paramter_recovery\param_recovery_test_data', 'param_recovery_test_data')


function [params,priors,rtbounds] = getAgentParams(agent)
%Previous function used to generate rts using optimal parameters

%This is frank specific but always return it to make life easy
priors.V=0; %initialize expected value for first trial to prior (possibly from previous run)
priors.Go=0; %initialize Go for first trial
priors.NoGo=0; %initialize NoGo for first trial
rtbounds = [1 5000]; %Don't let the agent choose 0

%Set agent parameters
switch agent
    case 'fixedLR_softmax'
        params = [0.03081894 .1 0.01698200]; %prop_spread  beta  alpha
    case 'fixedLR_egreedy'
        params = [0.45744891 0.15440092 0.11702531];
    case 'kalman_softmax'
        params=[0.5016875 0.53796360]; %Prop_spread Beta
    case 'kalman_processnoise'
        %         params=[0.44260383 0.2414416 2.384186e-07];
        params=[0.44260383 0.2414416 10]; %Changed omega to 10 to see if it could recover a non-zero number
    case 'kalman_sigmavolatility'
        %params = [0.272621000 0.10000000 0.00000000 0.3532713] prop_spread beta phi gamma
        params = [0.272621000 0.10000000 1 0.3532713]; %Changed phi to 1 to see if it could recover a non-zero number
        %params = [0.272621000 0.10000000 0 0]; %All condition params via recovery
    case 'franktc'
        params = [0.03059226 87511.530 3.449525 2.092848  685.891054 1879.997  611.3465]; %All condition params via Michael's email [lambda   epsilon   alphaG   alphaN  K   nu  rho]
    case 'kalman_uv_sum'
        %All condition params via Michael's email
        params =[0.62416337 0.001524063 0.69507251]; %prop_spread  beta  tau
    case 'kalman_uv_sum_kl'
        %All condition params via Michael's email
        params =[0.27298721 1.2574299 0.611199892 2.15673530 2.2536939758]; %prop_spread  beta  tau kappa lambda
    case 'kalman_uv_logistic'
        %All condition params via Michael's email
        params =[0.06442844 0.1209529 37.93337]; %prop_spread  tradeoff disrim
    case 'qlearning'
        %All condition params via Michael's email
        params =[0.98763,0.216741,0.1448303,0.9854957]; %gamma, alpha, epsilon lambda
    case 'sarsa'
        %All condition params via Michael's email
        params =[0.989816,0.239363,0.260227,0.9453975]; %gamma, alpha, epsilon lambda
    otherwise
        error('Not any agent I''ve heard of');
end


function [params,priors,rtbounds] = getAgentParamsUniform(agent,beta,n,a)
%Produce params via a uniform distribution each agent will be drawn from a
%specific seed


a_index = find(cellfun(@(x) strcmp(x,agent), {a.name}));
fprintf('%s\n',agent)
fprintf('a_index is %d\n',a_index)

%This is frank specific but always return it to make life easy
priors.V=0;                     %initialize expected value for first trial to prior (possibly from previous run)
priors.Go=0;                    %initialize Go for first trial
priors.NoGo=0;                  %initialize NoGo for first trial
rtbounds = [1 5000];            %Don't let the agent choose 0

%Fix the beta and some qlearning stats
beta = repmat(beta,n,1);
fix_gamma = .99;
fix_gamma = repmat(fix_gamma,n,1);
fix_lambda = .99;
fix_lambda = repmat(fix_lambda,n,1);

%Set agent parameters
switch agent
    case 'fixedLR_softmax'
        rng(1);
        prop_spread = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        alpha = genParamInRange(a(a_index).lower_bounds(2),a(a_index).upper_bounds(2),n);
        params = [prop_spread beta alpha]; %prop_spread  beta  alpha
    case 'fixedLR_egreedy' %Model not useful
        rng(2);
        params = []; 
    case 'kalman_softmax'
        rng(3);
        prop_spread = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        params=[prop_spread beta]; %Prop_spread Beta
    case 'kalman_processnoise'
        rng(4);
        prop_spread = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        omega = genParamInRange(a(a_index).lower_bounds(2),a(a_index).upper_bounds(2),n);
        params=[prop_spread beta omega]; %prop_spread beta omega
    case 'kalman_sigmavolatility' %Model not useful
        rng(5);
        params = []; %Changed phi to 1 to see if it could recover a non-zero number        
    case 'franktc'
        rng(6);
        lambda = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        epsilon = genParamInRange(a(a_index).lower_bounds(2),a(a_index).upper_bounds(2),n);
        alphaG = genParamInRange(a(a_index).lower_bounds(3),a(a_index).upper_bounds(3),n);
        alphaN = genParamInRange(a(a_index).lower_bounds(4),a(a_index).upper_bounds(4),n);
        K = genParamInRange(a(a_index).lower_bounds(5),a(a_index).upper_bounds(5),n);   
        nu = genParamInRange(a(a_index).lower_bounds(6),a(a_index).upper_bounds(6),n);  
        rho = genParamInRange(a(a_index).lower_bounds(7),a(a_index).upper_bounds(7),n);
        params = [lambda epsilon alphaG alphaN  K nu  rho]; %lambda   epsilon   alphaG   alphaN  K   nu  rho
    case 'kalman_uv_sum'
        rng(7);
        prop_spread = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        tau = genParamInRange(a(a_index).lower_bounds(2),a(a_index).upper_bounds(2),n);
        params =[prop_spread beta tau]; %prop_spread  beta  tau
    case 'kalman_uv_sum_kl' %Model not useful
        rng(8);
        params =[]; %prop_spread  beta  tau kappa lambda
    case 'kalman_uv_logistic'
        rng(9);
        prop_spread = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        tradeoff = genParamInRange(a(a_index).lower_bounds(2),a(a_index).upper_bounds(2),n);
        discrim = genParamInRange(a(a_index).lower_bounds(3),a(a_index).upper_bounds(3),n);
        params =[prop_spread tradeoff discrim]; %prop_spread  tradeoff disrim
    case 'qlearning'
        %Note that gamma and lambda are fixed, if this changes in the
        %initialize struct function the hard coded index numbers will have
        %to change
        rng(10);
        alpha = genParamInRange(a(a_index).lower_bounds(1),a(a_index).upper_bounds(1),n);
        epsilon = genParamInRange(a(a_index).lower_bounds(2),a(a_index).upper_bounds(2),n);
        params =[fix_gamma,alpha,epsilon,fix_lambda]; %gamma, alpha, epsilon lambda
    case 'sarsa' %Model not useful
        rng(11);
        params =[]; %gamma, alpha, epsilon lambda
    otherwise
        error('Not any agent I''ve heard of');

end

function param=genParamInRange(min,max,n)
%Generate random numbers within a specific range
%In:
%   min -> min range
%   max -> max range
%   n -> total number to generate
param = (max-min).*rand(n,1) + min;


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