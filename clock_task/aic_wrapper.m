agents = {'kalmanUV' 'kalmanLogistic' 'kalmanGRW' 'qlearning' 'sarsa'};
conditions = {'IEV' 'DEV' 'QUADUP' 'reversal'};
load('s.mat')
  for i = 1:length(agents)
      for j = 1:length(conditions)-1
          [s.(agents{i}).aic_stat(j,:),s.(agents{i}).aic_stat_trialwise(j,:)] = compute_aic(agents{i}, conditions{j},1,0);
      end
  end
  
  for o = 1:2
      rev_idx=o; %set  as 1 or 2,ItoD or DtoI
      iac_idx = rev_idx+3;
      for i = 1:length(agents)
          for j = 4
              [s.(agents{i}).aic_stat(iac_idx,:),s.(agents{i}).aic_stat_trialwise(iac_idx,:)] = compute_aic(agents{i}, conditions{j},1,rev_idx);
          end
      end
  end