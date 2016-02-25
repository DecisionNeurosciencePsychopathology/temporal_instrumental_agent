%Run this function to get error and bias infor on model parameter recovery
%Args are data,agent,and plot bool

function [T,models] = analyze_model_stability(fitted_vars,plots_on)

if nargin<2
    plots_on=0;
end

%Load in data
%load('fitted_vars')

%agent = 'kalman_uv_logistic'; %change angent name here


% org_params=getOrginalParams(agent);
% f_params = fitted_vars.noise_params.(agent).best_parameters;
% len = length(fitted_vars.noise_params.(agent).best_parameters);
% bias = mean(f_params) - org_params;
% bias_percent = bias./org_params.*100;
% SE = std(f_params)./sqrt(len);
% SE_percent = SE./org_params.*100;
% estimation_error=(sqrt(sum((repmat(org_params,len,1)-fitted_vars.noise_params.(agent).best_parameters).^2)))./(len);


%Tomorrow set this up so it will display everything???
%Fix the above code so it hits every model and save all that data in a
%table then use the plots on bool to determine if the histograms should be
%created  or not, essentially remove agent from args, loop over everything
%below, move plot on bool to above j loop...

models=fieldnames(fitted_vars.noise_params);
start_index= 1;
for i=1:length(models);
    model=models{i};
    num_params = size(fitted_vars.noise_params.(model).best_parameters,2);
    [org_params,param_names]=getOrginalParams(model);
    
    %param_names = tmp_param_names;
    
    
    %Calculate some descriptive stats
    f_params = fitted_vars.noise_params.(model).best_parameters;
    f_params_mean=mean(f_params);
    f_params_std=std(f_params);
    len = length(fitted_vars.noise_params.(model).best_parameters);
    bias = mean(f_params) - org_params;
    bias_percent = bias./org_params.*100;
    SE = std(f_params)./sqrt(len);
    SE_percent = SE./org_params.*100;
    estimation_error=(sqrt(sum((repmat(org_params,len,1)-fitted_vars.noise_params.(model).best_parameters).^2)))./(len);
    
    %Save in table for later
    T(start_index:start_index+num_params-1,:) = table(org_params',f_params_mean',f_params_std',bias',bias_percent',SE',SE_percent',estimation_error',...
        'RowNames',param_names,'VariableNames',{'original_parameters' 'fit_param_mean' 'fit_param_std' 'bias' 'bias_percent' 'SE' 'SE_percent' 'estimation_error'});
    
    %T.Properties.RowNames(start_index:start_index+num_params-1,1)=param_names;
    
    start_index = height(T) +1;
    
    
    
    param_str=[];
    if plots_on
        for j=1:num_params
            %fig_num = fig_num+1;
            figure(i)
            subplot(num_params,1,j)
            hist(fitted_vars.noise_params.(model).best_parameters(:,j),20);
            str=sprintf('Paramter space for %s %s org val: %.3f', model,param_names{j}, org_params(j));
            %str = ['Paramter space for ', model, ' ',param_names{j}, ' orginal parameter val: ' num2str(org_params(j))];
            title(str)
            param_str =[param_str '%.2f '];
        end
%         fprintf('Model: %s\n', model)
%         fprintf('Number of Parameters: %d\n', num_params)
%         fprintf(['Mean of Parameters: ',param_str,'\n\n'],mean(fitted_vars.noise_params.(model).best_parameters))
    end
end
end



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
end








% %Look at params via hist
% figure(88)
% hist(fitted_vars.noise_params.(agent).best_parameters(:,1));
% figure(89)
% hist(fitted_vars.noise_params.(agent).best_parameters(:,2));
% figure(90)
% hist(fitted_vars.noise_params.(agent).best_parameters(:,3));
%
%
% for j = 1:size(fitted_vars.noise_params.(agent).best_parameters,2)
%     figure(10+j)
%     hist(fitted_vars.noise_params.(agent).best_parameters(:,j),20);
% end