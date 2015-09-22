%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    behavfiles = glob('C:\kod\temporal_instrumental_agent\clock_task\subjects\*.csv');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'Alex')==1
        behavfiles = glob('/Users/localadmin/clock_smoothoperator/clock_task/subjects/*.csv');
    else
        behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
    end
end
sub=1;
behav{sub}.data = readtable(behavfiles{sub},'Delimiter',',','ReadVariableNames',true);
% write id
fname = behavfiles{sub};
idchars = regexp(fname,'\d');
behav{sub}.id = fname(idchars);
runs=unique(behav{sub}.data.run);

%load('fitted_vars_9-17-15.mat');
load('fitted_vars_noiseDiscount_alphaOmega.mat')
modelnames = fieldnames(fitted_vars);

for sub = 1:length(fitted_vars.fixed.best_parameters)
    for model = 1:numel(modelnames);
        
        modelname = modelnames{model};
        num_vars = length(fitted_vars.(modelname).best_parameters);
        %[aic_fixed(sub,model),bic_fixed(sub,model)] = aicbic(-fitted_vars.(modelname).best_cost(sub), num_vars, 400);
        [aic(sub,model),bic(sub,model)] = aicbic(-fitted_vars.(modelname).best_cost(sub), num_vars, 400);
        
    end
    
    
    
end


for m = 1:model
    modelname = modelnames{m};
    %Save sum bic
    fitted_vars.(modelname).sum_bic = sum(bic(:,m));
    
    
    %Save sum costs
    fitted_vars.(modelname).sum_cost = sum(fitted_vars.(modelname).best_cost);
    temp_costs(:,m) =  sum(fitted_vars.(modelname).best_cost);
    
    
    %Plot the data for easy viewing
    figure(88);
    hold on
    plot(1:sub, aic(:,m))
    title('AIC comparision of all models')
    legend(modelnames{:});
    
    figure(89);
    hold on
    plot(1:sub, bic(:,m))
    title('BIC comparision of all models')
    legend(modelnames{:});
    
    figure(90);
    hold on
    plot(1:sub, fitted_vars.(modelname).best_cost)
    title('Comparision of costs between all models')
    legend(modelnames{:});
end

%Print some statistics
[val_BIC,lowest_BIC] = min(sum(bic));
[val_costs,lowest_cost] = min(temp_costs);
fprintf('Lowest overall BIC: %d model: %s\n',val_BIC, modelnames{lowest_BIC});
fprintf('Lowest overall cost: %d model: %s\n',val_costs, modelnames{lowest_cost});
