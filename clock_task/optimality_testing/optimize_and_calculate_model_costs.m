%Script to optimize parameters of RL models and calculate costs over a
%specified number of runs

%%
%basic structure
%cd '/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/optimality_testing';
agents = initialize_agents_struct;
%agentnames = fieldnames(agents);
%nagents = length(agentnames);
nagents = length(agents);
noptim = 100; %optimize parameters 100x
ncosts = 1000;
nruns=10; %10 runs per contingency
ntrials = 150;
nbasis = 24;
ntimesteps = 500;
load('mastersamp.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
condnames=fieldnames(mastersamp); %by default, optimize over all contingencies in lookup
ncond=length(condnames); %how many conditions/contingencies

addpath('../');

%truncate mastersamp to the number of trials used here to save RAM
for i = 1:length(condnames)
    mastersamp.(condnames{i}).lookup = mastersamp.(condnames{i}).lookup(:,1:ntrials);
end

%generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
%handle permutation outside of optimization for speed
optmat = cell(length(condnames),1);
for i = 1:ncond
    clear row; %need variable to be deleted for array of structs below to work
    for j = 1:nruns
        tmp = mastersamp.(condnames{i});
        tmp.name = condnames{i}; %for keeping track of what contingency this is (esp. when they mix in ALL)
        tmp.lookup = tmp.lookup(:, randperm(size(tmp.lookup,2))); %randomly permute columns
        %tmp.perm = randperm(size(tmp.lookup,2)); %just for checking that this works
        %tmp.lookup = tmp.lookup(:, tmp.perm); %randomly permute columns
        row(j) = tmp; %add to a vector of structs
    end
    optmat{i} = row;
end

%concatenate optmat to obtain an ALL condition
allcond = horzcat(optmat{:});
optmat{length(optmat) + 1} = allcond;
condnames{length(condnames) + 1} = 'ALL';
ncond = ncond + 1; %add all condition

%optmat is now a cell vector where each element is an vector of structs for the lookups to be used for each run.
%the final element contains all of the preceding conditions. For the five contingencies and 10 runs, this amounts
%to a 58MB object that would need to be passed along to parpool workers.

%add opt_results and opt_costs fields to each agent (based on number of optimizations)
for i = 1:nagents
    agents(i).nbasis = nbasis;
    agents(i).ntimesteps = ntimesteps;
    agents(i).ntrials = ntrials;
    agents(i).runseed = 888; %seeds rngs inside runs?
    %agents(i).opt_params = NaN(noptim, ncond, agents(i).k);
    %agents(i).opt_costs = NaN(noptim, ncond);
    %agents(i).opt_conds = condnames;
    
    %agents.(agentnames{i}).opt_params = NaN(noptim, ncond, agents.(agentnames{i}).k);
    %agents.(agentnames{i}).opt_costs = NaN(noptim, ncond);
    %agents.(agentnames{i}).opt_conds = condnames;
end

clear mastersamp tmp row allcond i j;

%so, MATLAB requires that for arrays defined outside of loop, but operated on within the parfor,
%the first-level index must be a variant of the parfor counter (i).
%this essentially kills the possibility of loop unrolling because i is never used directly as a first-level index!

costs=NaN(noptim, nagents, ncond);
pars=cell(noptim, nagents, ncond);
poolobj=parpool('local',40);
parfor i = 1:noptim
%for i = 1:noptim
    ipars=cell(nagents, ncond);
    icosts=NaN(nagents, ncond);
    
    for a = 1:nagents
        for c = 1:ncond
            [ipars{a, c}, icosts(a, c)] = GAoptimize(agents(a), optmat{c});
        end
    end
   
    costs(i, :, :) = icosts;
    pars(i, :, :) = ipars;
end
delete(poolobj);
save('optimize_output.mat', 'costs', 'pars', 'agents', 'optmat');

%assign pars and costs back into agents array?

%%beautiful, but impossible in matlab. leaving here to cry over

%optimize parameters in parallel.
%need to unroll loop over optimization runs and agents for maximal speed.
%matlab does not allow nested parfor loops, so we need parfor to span out to noptim x nagents.
% nloops = noptim*nagents;
% 
% parfor i = 1:nloops
%     %divide into blocks of length noptim (100)
%     if mod(i-1, noptim) == 0
%         a = (i-1)/noptim + 1;
%     end
%     
%     optnum = mod(i-1, noptim) + 1; %1..100
%     %optresults = NaN(ncond,agents(a).k);
%     optcosts = NaN(ncond,1);
%     
%     %loop over contingencies, optimizing parameters
%     for j = 1:ncond
%         [optresults(j,:), optcosts(j)] = GAoptimize(agents(a), optmat{j});
%     end
%     
%     agents(a).opt_params(optnum,:,:) = optresults; %assign opimal parameters for all conditions of the current optimization to struct
%     agents(a).opt_costs(optnum,:) = optcosts; 
% end

    




% %Load in main struct s
% load('s.mat');
% 
% %Load in clock options
% load('clock_options.mat');
% 
% %load in lookup tables
% load('mDEV.mat');
% load('mIEV.mat');
% load('mQUADUP.mat');
% 
% lookups = {mIEV mDEV mQUADUP}; %replaces cond
% agents = fieldnames(s);
% optimize=0; %set to calc contingency indivisually
% optimizeAll=1; %set to 1 to calc 'all' contingencies
% calc_cost=0; %Set to 1 to calulate costs of nruns
% calc_cost_all=0; %1 = find cost using 'all condition' optimal params
% calc_cost_reversal=0; %0 = off; 1 = IEV to DEV; 2 = DEV to IEV
% 
% %Trial set up variables
% ntrials = 200;
% nbasis = 24;
% ntimesteps = 500;
% 
% 
% %Set up seeds
% GAnruns = 15; %This is the nruns variable in the fitness function
% nruns = 100;
% runseed = 55;
% rng(runseed);
% seeds=randi([1 500], nruns, 4);
% 
% 
% %Set up permutation of reward martix here
% global vperm
% vperm = zeros(GAnruns,length(mIEV.sample));
% for i = 1:GAnruns
%     vperm(i,:) = randperm(length(mIEV.lookup));
% end
% 
% for i=1:length(lookups)+1
%     
%     if i<=3
%         m = lookups{i};
%         mAll{1,i} = m;
%     end
%     if optimize
%         for j = 1:length(agents)
%             %Use set field here to capture opt values
%             [s.(agents{j}).opt_params(i,:), s.(agents{j}).Fga(i,:)] = GAoptimize(s.(agents{j}).init_params,...
%                 s.(agents{j}).lower_bounds,s.(agents{j}).upper_bounds,s.(agents{j}).numVars,j,m);
%         end
%     elseif optimizeAll && i>3
%         for j = 1:length(agents)
%             %Use set field here to capture opt values
%             [s.(agents{j}).opt_params(i,:), s.(agents{j}).Fga(i,:)] = GAoptimize(s.(agents{j}).init_params,...
%                 s.(agents{j}).lower_bounds,s.(agents{j}).upper_bounds,s.(agents{j}).numVars,j,mAll);
%         end
%         
%         %Save param data in temp struct...Power outages :(
%         s_tmp = s;
%         save s_tmp s_tmp;
%         
%         if calc_cost_all
%             break; %Kick out for now until I get all optimal params again
%         end
%         
%     end
%     
%     %Calculate the condition specific costs using condition specific
%     %params
%     if(calc_cost) && i<=3 %indexing issue
%         for o = 1:nruns
%             rev=0;
%             vperm_run = randperm(length(mIEV.sample));
%             m.lookup = m.lookup(:,vperm_run); %permute here as well
%             fprintf('Computing costs with condition specific optimal params, run number %d and rngseeds: %s \n', o, num2str(seeds(o,:)));
%             s = calcCost(s,agents,seeds,i,m,o,ntrials,nbasis,ntimesteps,i,clock_options,'costs',rev);
%         end
%     end
%     
%     %Calculate the costs for the 'all conditions' optimal parameters
%     if(calc_cost_all) && i<=3 %indexing issue
%         for o = 1:nruns
%             rev=0;
%             vperm_run = randperm(length(mIEV.sample));
%             m.lookup = m.lookup(:,vperm_run); %permute here as well
%             fprintf('Computing costs with All optimal params, run number %d and rngseeds: %s \n', o, num2str(seeds(o,:)));
%             s = calcCost(s,agents,seeds,i,m,o,ntrials,nbasis,ntimesteps,4,clock_options,'costsAll',rev);
%         end
%     end
%     
%     %Using the condition all optimal parmas calculate the cost with a
%     %reversal, IEV to DEV
%     if(calc_cost_reversal==1) && i==1
%         for o = 1:nruns
%             rev=1;
%             vperm_run = randperm(length(mIEV.sample));
%             m.vperm_run = vperm_run;
%             m.lookup = m.lookup(:,vperm_run); %permute here as well
%             fprintf('Computing costs with All optimal params, and reversal, run number %d and rngseeds: %s \n', o, num2str(seeds(o,:)));
%             s = calcCost(s,agents,seeds,i,m,o,ntrials,nbasis,ntimesteps,4,clock_options,'costsIevtoDev', rev);
%             %s.kalmanSKEPTIC.costsIevtoDev(i,o) %debug check
%         end
%     end
%     
%     %DEV to IEV
%     if(calc_cost_reversal==2) && i==2
%         for o = 1:nruns
%             rev=1;
%             vperm_run = randperm(length(mIEV.sample));
%             m.vperm_run = vperm_run;
%             m.lookup = m.lookup(:,vperm_run); %permute here as well
%             fprintf('Computing costs with All optimal params, and reversal, run number %d and rngseeds: %s \n', o, num2str(seeds(o,:)));
%             s = calcCost(s,agents,seeds,i,m,o,ntrials,nbasis,ntimesteps,4,clock_options,'costsDevtoIev', rev);
%             %s.kalmanSKEPTIC.costsIevtoDev(i,o) %debug check
%         end
%     end
%     
% end
% 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is just used for initialization or running some simple
%stats on the structure produced from above

%Used this to set up s
% agents = {'kalmanSKEPTIC' 'kalmanLogistic' 'kalmanGRW' 'qlearning' 'sarsa'};
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
%  agents = {'kalmanSKEPTIC' 'kalmanLogistic' 'kalmanGRW' 'qlearning' 'sarsa'};
%  for i = 4:length(agents)
%      for w = 1:3 %For each condition
% % % %         s.(agents{i}).mean_cost(w,:)=mean(s.(agents{i}).costs(w,:),2);
% % % %         s.(agents{i}).sum_cost(w,:)=sum(s.(agents{i}).costs(w,:),2);
% % % %         s.(agents{i}).std_cost(w,:)=std(s.(agents{i}).costs(w,:),0,2);
% % % %         s.(agents{i}).max_cost(w,:)=max(s.(agents{i}).costs(w,:));
% % % %         s.(agents{i}).min_cost(w,:)=min(s.(agents{i}).costs(w,:));
% % % %
% % %          %Compute costs stats using 'All' optimal params
%          s.(agents{i}).mean_costAll(w,:)=mean(s.(agents{i}).costsAll(w,:),2);
%          s.(agents{i}).sum_costAll(w,:)=sum(s.(agents{i}).costsAll(w,:),2);
%          s.(agents{i}).std_costAll(w,:)=std(s.(agents{i}).costsAll(w,:),0,2);
%          s.(agents{i}).max_costAll(w,:)=max(s.(agents{i}).costsAll(w,:));
%          s.(agents{i}).min_costAll(w,:)=min(s.(agents{i}).costsAll(w,:));
%
%
%          s.(agents{i}).mean_costIevtoDev(1,:)=mean(s.(agents{i}).costsIevtoDev(1,:));
%
%
%          s.(agents{i}).mean_costDevtoIev(1,:)=mean(s.(agents{i}).costsDevtoIev(2,:));
%
%
%
%      end
% end

