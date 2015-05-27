function iac_stat=compute_iac(model, condition, all)
%This script is to calculate the iac statistic for the clock models (kalman, TD, frank)
%n is the number of trials
%k is the number of params

load('s.mat');
load('perfect_cost.mat');
n=200; %number of trials


if (all) val = 'mean_costAll'; else val = 'mean_cost'; end

switch condition
    case 'IEV'
        idx = 1;
    case 'DEV'
        idx = 2;
    otherwise
        idx = 3;
end

if strfind(model, 'kalman') %will have to update once Frank is introduced
    k=2;
else
    k=4;
end


model_cost = s.(model).(val)(idx);
 
%The equation
iac_stat = log(((perfect_cost.(condition) - model_cost).^2)/n)+2*k;