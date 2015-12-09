%Remove outliers or bad subjects
%Note: Do make a backup just incase you delete the wrong person

load('param_recov.mat')
agents = fieldnames(fitted_vars.subj_fitting);
%Subject to exclude index numer 1-77 currently
exclude_idx = 24;
for j = 1:length(exclude_idx)
    for i = 1:length(agents)
        fitted_vars.subj_fitting.(agents{i}).best_cost(exclude_idx(j))=[];
        fitted_vars.subj_fitting.(agents{i}).best_parameters(exclude_idx(j),:)=[];
    end
end


save param_recov fitted_vars