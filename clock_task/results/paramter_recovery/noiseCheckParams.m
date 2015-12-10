% std_mean_error=std(fitted_vars.noise_params.fixed_rho.best_parameters)/sqrt(100);

agent = 'kalman_uv_logistic'; %change angent name here

org_params = [0.62416337 0.001524063 0.69507251];
f_params = fitted_vars.noise_params.(agent).best_parameters;
len = length(fitted_vars.noise_params.(agent).best_parameters);
bias = mean(f_params) - org_params;
bias_percent = bias./org_params.*100;
SE = std(f_params)./sqrt(len);
SE_percent = SE./org_params.*100;
estimation_error=(sqrt(sum((repmat(org_params,len,1)-fitted_vars.noise_params.(agent).best_parameters).^2)))./(len);
%sqrt(sum((0.2-fitted_vars.noise_params.(agent).best_parameters(:,1)).^2)/100)







models=fieldnames(fitted_vars.noise_params); 
for i=1:length(models);
    model=models{i};
    num_params = size(fitted_vars.noise_params.(model).best_parameters,2);
    param_str=[];
    for j=1:num_params
       %fig_num = fig_num+1; 
       figure(i)
       subplot(num_params,1,j)
       hist(fitted_vars.noise_params.(model).best_parameters(:,j),20);
       title(['Paramter space for ', model, ' parameter', num2str(j)])
       param_str =[param_str '%.2f ']; 
    end
    fprintf('Model: %s\n', model)
    fprintf('Number of Parameters: %d\n', num_params)
    fprintf(['Mean of Parameters: ',param_str,'\n\n'],mean(fitted_vars.noise_params.(model).best_parameters))
end




% function org_params=getOrginalParams(agent)
% switch agent
%     case 'fixedLR_softmax'
%         params = [0.03081894 .1 0.01698200]; %prop_spread  beta  alpha
%     case 'fixedLR_egreedy'
%         params = [0.45744891 0.15440092 0.11702531];
%     case 'kalman_softmax'
%         params=[0.5016875 0.53796360]; %Prop_spread Beta
%     case 'kalman_processnoise'
%         params=[0.44260383 0.2414416 2.384186e-07];
%     case 'kalman_sigmavolatility'
%         params = [0.272621000 0.10000000 0.00000000 0.3532713]; %All condition params via Michael's email
%         %params = [0.272621000 0.10000000 0 0]; %All condition params via recovery
%     case 'frank'
%         params = [0.03059226 87511.530 3.449525 2.092848  685.891054 1879.997  611.3465]; %All condition params via Michael's email [lambda   epsilon   alphaG   alphaN  K   nu  rho]
%         priors.V=0; %initialize expected value for first trial to prior (possibly from previous run)
%         priors.Go=0; %initialize Go for first trial
%         priors.NoGo=0; %initialize NoGo for first trial
%         rtbounds = [1 5000]; %Don't let the agent choose 0
%     case 'kalman_uv_sum'
%         %All condition params via Michael's email
%         params =[0.62416337 0.001524063 0.69507251]; %prop_spread  beta  discrim
%     case 'kalman_uv_sum_kl'
%         %All condition params via Michael's email
%         params =[0.27298721 1.2574299 0.611199892 2.15673530 2.2536939758]; %prop_spread  beta  tau kappa lambda
%     case 'kalman_uv_logistic'
%         %All condition params via Michael's email
%         params =[0.62416337 0.001524063 0.69507251]; %prop_spread  tradeoff disrim
%     case 'qlearning'
%         %All condition params via Michael's email
%         params =[0.98763,0.216741,0.1448303,0.9854957]; %gamma, alpha, epsilon lambda
%     case 'sarsa'
%         %All condition params via Michael's email
%         params =[0.989816,0.239363,0.260227,0.9453975]; %gamma, alpha, epsilon lambda
%     otherwise
%         error('Not any agent I''ve heard of');
% end
% end








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