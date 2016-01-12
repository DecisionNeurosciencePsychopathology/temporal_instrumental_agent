%Compute the AIC and BIC of the fixed and fixed discounted (fixed_KL)
%models


%load('fitted_vars_9-17-15.mat');
%load('fitted_vars_noise_run_1.mat');
load('param_recov.mat')
%fitted_vars=rmfield(fitted_vars,'noise_params');
%fitted_vars=rmfield(fitted_vars,'optimal_params');
fit_type = fieldnames(fitted_vars);
modelnames = fieldnames(fitted_vars.subj_fitting);
pick = [1 2 3 4 5 6 7]; % skip Q and Frank for now;
modelnames = modelnames(pick);
%modelnames = {'fixed' 'fixed_KL'};
%num_vars = [1 3];
%load('fitted_vars_9-16-15.mat');
n = 400;

for i = 1:length(fit_type)
    fit = fit_type{i};
    for sub = 1:length(fitted_vars.(fit).fixedLR_softmax.best_cost)
        for model = 1:numel(modelnames);
            
            modelname = modelnames{model};
            num_vars = size(fitted_vars.subj_fitting.(modelname).best_parameters,2);
            pseudoLogL = -n.*log(sqrt(fitted_vars.(fit).(modelname).best_cost(sub))./n)./2;  
            fitted_vars.subj_fitting.(modelname).aic = pseudoLogL;
            %[aic_fixed(sub,model),bic_fixed(sub,model)] = aicbic(-fitted_vars.(modelname).best_cost(sub), num_vars, 400);
%             [aic(sub,model),bic(sub,model)] = aicbic(pseudoLogL, num_vars, n);
            [aic(sub,model)] = aicbic(pseudoLogL, num_vars, n);

        end
    end
    [p,tbl,stats] = anova1(aic,modelnames);
    title('AIC Boxplot');

    
    %Plot the data for easy viewing
    figure(90+i); clf;
    plot(1:sub, aic)
    hold on
    %plot(1:sub, aic(:,2))
    %plot(1:sub, aic(:,3))
    title('AIC comparision by model and subject')
    %legend('fixed','fixedKL', 'noiseDiscount');
    legend(modelnames);
    
%     figure(100+i); clf;
%     plot(1:sub, bic)
%     hold on
%     % plot(1:sub, bic(:,2), 'r')
%     % plot(1:sub, aic(:,3), 'g')
%     title('BIC comparision of fixed and fixed KL')
%     %legend('fixed','fixedKL', 'noiseDiscount');
%     legend(modelnames);
    
    %Use smp_BMS to comput the Bayesian model slection for group studies  Bayesian model selection for group studies
    %   FORMAT [alpha, exp_r, xp] = spm_BMS (lme, Nsamp, do_plot, sampling, ecp, alpha0)
    %
    %   INPUT:
    %   lme      - array of log model evidences
    %                rows: subjects
    %                columns: models (1..Nk)
    %   Nsamp    - number of samples used to compute exceedance probabilities
    %              (default: 1e6)
    %   do_plot  - 1 to plot p(r|y)
    %   sampling - use sampling to compute exact alpha
    %   ecp      - 1 to compute exceedance probability
    %   alpha0   - [1 x Nk] vector of prior model counts
    %
    %   OUTPUT:
    %   alpha   - vector of model probabilities
    %   exp_r   - expectation of the posterior p(r|y)
    %   xp      - exceedance probabilities
    [alpha, exp_r, xp ] = spm_BMS(-aic,1e6)
    figure(200+i); clf;
    subplot(2,1,1);
    bar(exp_r); 
    set(gca, 'XTick', 1:length(modelnames), 'XTickLabel', modelnames);
    title('Expectation of the posterior')
    subplot(2,1,2);
    bar(xp);
    set(gca, 'XTick', 1:length(modelnames), 'XTickLabel', modelnames);
    title('Exceedance probability')
end