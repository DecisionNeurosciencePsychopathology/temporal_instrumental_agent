%Script to optimize parameterds of RL models and calculate costs over a
%given amount of runs

%Load in main struct s
load('s.mat');

%Load in clock options
load('clock_options.mat');

%load in lookup tables
load('mDEV.mat');
load('mIEV.mat');
load('mQUADUP.mat');

lookups = {mIEV mDEV mQUADUP}; %replaces cond
agents = fieldnames(s);
optimize=0; %set to calc contingency indivisually
optimizeAll=1; %set to 1 to calc 'all' contingencies
calc_cost=0; %Set to 1 to calulate costs of nruns
calc_cost_all=0;

%Trial set up variables
ntrials = 200;
nbasis = 24;
ntimesteps = 500;


%Set up seeds
GAnruns = 15; %This is the nruns variable in the fitness function
nruns = 100;
runseed = 55;
rng(runseed);
seeds=randi([1 500], nruns, 4);

%Set up permutation of reward martic here
global vperm 
vperm = zeros(GAnruns,length(mIEV.sample));
for i = 1:GAnruns
    vperm(i,:) = randperm(length(mIEV.lookup));
end

for i=1:length(lookups)+1
    if i<=3
        m = lookups{i};
        mAll{1,i} = m;
    end
    if optimize
        for j = 1:length(agents)
            %Use set field here to capture opt values
            [s.(agents{j}).opt_params(i,:), s.(agents{j}).Fga(i,:)] = GAoptimize(s.(agents{j}).init_params,...
                s.(agents{j}).lower_bounds,s.(agents{j}).upper_bounds,s.(agents{j}).numVars,j,m);
        end
    elseif optimizeAll && i>3
        for j = 1:length(agents)
            %Use set field here to capture opt values
            [s.(agents{j}).opt_params(i,:), s.(agents{j}).Fga(i,:)] = GAoptimize(s.(agents{j}).init_params,...
                s.(agents{j}).lower_bounds,s.(agents{j}).upper_bounds,s.(agents{j}).numVars,j,mAll);
        end
        if calc_cost_all
            break; %Kick out for now until I get all optimal params again
        end
    end
    
        %This loop will compute the costs of the agents
        if(calc_cost) && i<=3 %indexing issue
            for o = 1:nruns
                vperm_run = randperm(length(mIEV.sample));
                m.lookup = m.lookup(:,vperm_run); %permute here as well
                fprintf('Computing costs with condition specific optimal params, run number %d and rngseeds: %s \n', o, num2str(seeds(o,:)));
                s.(agents{1}).costs(i,o) = clock_logistic_operator_kalman_stats(s.(agents{1}).opt_params(i,:),seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 1, m);
                s.(agents{2}).costs(i,o) = clock_logistic_operator_kalman_stats([s.(agents{2}).opt_params(i,:) 0 0.2],seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 0, m);
                s.(agents{3}).costs(i,o) = clock_logistic_operator_kalman_stats([s.(agents{3}).opt_params(i,:) 0.9 0.2],seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 0, m);
                
                %Specify agent
                clock_options.agent = 'qlearning';
                s.(agents{4}).costs(i,o) = ClockWalking_3D_discountedEv_stats(clock_options,m, seeds(o,1:2),s.(agents{4}).opt_params(i,:));
                clock_options.agent = 'sarsa';
                s.(agents{5}).costs(i,o) = ClockWalking_3D_discountedEv_stats(clock_options,m, seeds(o,1:2),s.(agents{5}).opt_params(i,:));
            end
        end
        
        if(calc_cost_all) && i<=3 %indexing issue 
            for o = 1:nruns
                vperm_run = randperm(length(mIEV.sample));
                m.lookup = m.lookup(:,vperm_run); %permute here as well
                fprintf('Computing costs with All optimal params, run number %d and rngseeds: %s \n', o, num2str(seeds(o,:)));
                s.(agents{1}).costsAll(i,o) = clock_logistic_operator_kalman_stats(s.(agents{1}).opt_params(4,:),seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 1, m);
                s.(agents{2}).costsAll(i,o) = clock_logistic_operator_kalman_stats([s.(agents{2}).opt_params(4,:) 0 0.2],seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 0, m);
                s.(agents{3}).costsAll(i,o) = clock_logistic_operator_kalman_stats([s.(agents{3}).opt_params(4,:) 0.9 0.2],seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 0, m);
                
                %Specify agent
                clock_options.agent = 'qlearning';
                s.(agents{4}).costsAll(i,o) = ClockWalking_3D_discountedEv_stats(clock_options,m, seeds(o,1:2),s.(agents{4}).opt_params(4,:));
                clock_options.agent = 'sarsa';
                s.(agents{5}).costsAll(i,o) = ClockWalking_3D_discountedEv_stats(clock_options,m, seeds(o,1:2),s.(agents{5}).opt_params(4,:));
            end
        end
end


% for p=1:5
%     s.(agents{p}).costsAll =0;
% end



%This code is just used for initialization or running some simple
%stats on the structure produced from above

%Used this to set up s
% agents = {'kalmanUV' 'kalmanLogistic' 'kalmanGRW' 'qlearning' 'sarsa'};
% for i = 1:length(agents)
%     s.(agents{i}).init_params=0;
%     s.(agents{i}).lower_bounds=0;
%     s.(agents{i}).upper_bounds=0;
%     s.(agents{i}).numVars=0;
%     s.(agents{i}).opt_params=0;
%     s.(agents{i}).Fga=0;
%     s.(agents{i}).costs=0;
%     s.(agents{i}).costsAll=0;
% end
% 
% for i = 1:3
%     s.(agents{i}).init_params=[0.6 0.08];
%     s.(agents{i}).lower_bounds=[0.1 0.01];
%     s.(agents{i}).upper_bounds=[0.99 0.25];
%     s.(agents{i}).numVars=2;
% s.(agents{i}).opt_params=zeros(1,s.(agents{i}).numVars);
% end
% 
% for i = 4:5
%     s.(agents{i}).init_params=[0.9 0.2 0.08 0.99];
%     s.(agents{i}).lower_bounds=[0.8 0.01 0.01 0.9];
%     s.(agents{i}).upper_bounds=[0.99 0.35 0.3 0.999];
%     s.(agents{i}).numVars=4;
% s.(agents{i}).opt_params=zeros(1,s.(agents{i}).numVars);
% end


%Used this to set up stats for models
% agents = {'kalmanUV' 'kalmanLogistic' 'kalmanGRW' 'qlearning' 'sarsa'};
% for i = 1:length(agents)
%     for w = 1:3 %For each condition
%         s.(agents{i}).mean_cost(w,:)=mean(s.(agents{i}).costs(w,:),2);
%         s.(agents{i}).sum_cost(w,:)=sum(s.(agents{i}).costs(w,:),2);
%         s.(agents{i}).std_cost(w,:)=std(s.(agents{i}).costs(w,:),0,2);
%         s.(agents{i}).max_cost(w,:)=max(s.(agents{i}).costs(w,:));
%         s.(agents{i}).min_cost(w,:)=min(s.(agents{i}).costs(w,:));
% 
%         %Compute costs stats using 'All' optimal params
%         s.(agents{i}).mean_costAll(w,:)=mean(s.(agents{i}).costsAll(w,:),2);
%         s.(agents{i}).sum_costAll(w,:)=sum(s.(agents{i}).costsAll(w,:),2);
%         s.(agents{i}).std_costAll(w,:)=std(s.(agents{i}).costsAll(w,:),0,2);
%         s.(agents{i}).max_costAll(w,:)=max(s.(agents{i}).costsAll(w,:));
%         s.(agents{i}).min_costAll(w,:)=min(s.(agents{i}).costsAll(w,:));
%     end
% end


