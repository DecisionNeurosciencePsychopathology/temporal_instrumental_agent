function [p,tbl,stats] =see_subject_params(fitted_vars)

models=fieldnames(fitted_vars.subj_fitting); 
pick = [1 2 3 6]; % skip Q and Frank for now;
models = models(pick);
y=[];
%Plotting loop to show subject parameter space for each model
for i=1:length(models);
    model=models{i};
    [~,names]=getOrginalParams(model);
    num_params = size(fitted_vars.subj_fitting.(model).best_parameters,2);
    param_str=[];
    for j=1:num_params
       %fig_num = fig_num+1; 
       figure(i)
       subplot(num_params,1,j)
       hist(fitted_vars.subj_fitting.(model).best_parameters(:,j),20);
       title(['Paramter space for ', model, ' parameter ', names{j}])
       param_str =[param_str '%.2f '];
    end
      if ~strcmp(model,'kalman_softmax')
       subplot(num_params+1,1,j)

       scatter(fitted_vars.subj_fitting.(model).best_parameters(:,j),fitted_vars.subj_fitting.(model).best_parameters(:,j-1));
       xlabel(['Parameter: ', names{j}]);
       ylabel(['Parameter: ', names{j-1}]);
      end
    fprintf('Model: %s\n', model)
    fprintf('Number of Parameters: %d\n', num_params)
    fprintf(['Mean of Parameters: ',param_str,'\n\n'],mean(fitted_vars.noise_params.(model).best_parameters))
    
    if strcmp(model,'qlearning')
        factr = 100;
    else
        factr=1;
    end
    
    y(:,i) = fitted_vars.subj_fitting.(model).best_cost.*factr;
end

[p,tbl,stats] = anova1(y,models);




function [params,names]=getOrginalParams(agent)
%JW 11/11:
%Beta's were fixed at .1 at the time of writing this, also gamma and lambda
%was fixed as well
switch agent
    case 'fixedLR_softmax'
        params = [0.03081894 0.01698200]; %prop_spread  beta  alpha
        names={'prop_spread'; 'alpha'};
    case 'fixedLR_egreedy'
        params = [0.45744891 0.15440092 0.11702531]; %Prop_spread epsilon alpha
        names={'prop_spread';'epslion';'alpha'};
    case 'kalman_softmax'
        params=[0.5016875]; %Prop_spread Beta
        names={'prop_spread'};
    case 'kalman_processnoise'
        %params=[0.44260383 2.384186e-07]; %Prop_spread Beta Omega
        params=[0.44260383  10];
        names={'prop_spread'; 'omega'};
    case 'kalman_sigmavolatility'
        %params = [0.272621000 0.00000000 0.3532713]; %prop_spread   phi     gamma
        params = [0.272621000 1 0.3532713];
        names={'prop_spread'; 'phi'; 'gamma'};
        
    case 'franktc'
        params = [0.03059226 87511.530 3.449525 2.092848  685.891054 1879.997  611.3465]; %All condition params via Michael's email [lambda   epsilon   alphaG   alphaN  K   nu  rho]
        names={'lambda';   'epsilon';   'alphaG';   'alphaN';  'K';   'nu';  'rho';};
    case 'kalman_uv_sum'
        %All condition params via Michael's email
        params =[0.62416337 0.69507251]; %prop_spread  beta  discrim
        names={'prop_spread';'discrim'};
    case 'kalman_uv_sum_kl'
        %All condition params via Michael's email
        params =[0.27298721 0.611199892 2.15673530 2.2536939758]; %prop_spread  beta  tau kappa lambda
        names={'prop_spread'; 'tau'; 'kappa'; 'lambda';};
    case 'kalman_uv_logistic'
        %All condition params via Michael's email
        %params =[0.62416337 0.001524063 0.69507251]; %prop_spread  tradeoff disrim
        params =[0.06442844 0.1209529 37.93337];
        names={'prop_spread'; 'tradeoff'; 'discrim'};
    case 'qlearning'
        %All condition params via Michael's email
        %params =[0.98763,0.216741,0.1448303,0.9854957]; %gamma, alpha, epsilon lambda
        params =[0.216741,0.1448303]; %gamma, alpha, epsilon lambda
        names={'alpha'; 'epsilon'};
    case 'sarsa'
        %All condition params via Michael's email
        %params =[0.989816,0.239363,0.260227,0.9453975]; %gamma, alpha, epsilon lambda
        params =[0.239363,0.260227]; %gamma, alpha, epsilon lambda
        names={'alpha'; 'epsilon';};
    otherwise
        error('Not any agent I''ve heard of');
end