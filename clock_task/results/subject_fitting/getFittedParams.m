function [opt_params,search_params,names]=getFittedParams(fitted_vars,sub,agent)
beta = .1;

switch agent
    case 'fixedLR_softmax'
        opt_params = [0.03081894 1.536380918 0.01698200]; %prop_spread  beta  alpha
        search_params = [fitted_vars.subj_fitting.(agent).best_parameters(sub,1) beta fitted_vars.subj_fitting.(agent).best_parameters(sub,2)];
        names={'prop_spread'; 'alpha'};
    case 'fixedLR_egreedy'
        opt_params = [0.45744891 0.15440092 0.11702531]; %Prop_spread epsilon alpha
        names={'prop_spread';'epslion';'alpha'};
    case 'kalman_softmax'
        opt_params=[0.5016875 0.53796360]; %Prop_spread Beta
        search_params = [fitted_vars.subj_fitting.(agent).best_parameters(sub,1) beta];
        names={'prop_spread'};
    case 'kalman_processnoise'
        %params=[0.44260383 2.384186e-07]; %Prop_spread Beta Omega
        opt_params=[0.44260383 0.2414416 2.384186e-07];
        %opt_params=[0.44260383 0.2414416 10];
        search_params = [fitted_vars.subj_fitting.(agent).best_parameters(sub,1) beta fitted_vars.subj_fitting.(agent).best_parameters(sub,2)];
        names={'prop_spread'; 'omega'};
    case 'kalman_sigmavolatility'
        %params = [0.272621000 0.00000000 0.3532713]; %prop_spread   phi     gamma
        opt_params = [0.272621000 1 0.3532713];
        names={'prop_spread'; 'phi'; 'gamma'};
    case 'franktc'
        opt_params = [0.03059226 87511.530 3.449525 2.092848  685.891054 1879.997  611.3465]; %All condition params via Michael's email [lambda   epsilon   alphaG   alphaN  K   nu  rho]
        search_params = fitted_vars.subj_fitting.franktc.best_parameters;
        names={'lambda';   'epsilon';   'alphaG';   'alphaN';  'K';   'nu';  'rho';};
    case 'kalman_uv_sum'
        %All condition params via Michael's email
        opt_params =[0.62416337 0.001524063 0.69507251]; %prop_spread  beta  discrim
        search_params = [fitted_vars.subj_fitting.(agent).best_parameters(sub,1) beta fitted_vars.subj_fitting.(agent).best_parameters(sub,2)];
        names={'prop_spread';'discrim'};
    case 'kalman_uv_sum_kl'
        %All condition params via Michael's email
        opt_params =[0.27298721 0.611199892 2.15673530 2.2536939758]; %prop_spread  beta  tau kappa lambda
        names={'prop_spread'; 'tau'; 'kappa'; 'lambda';};
    case 'kalman_uv_logistic'
        %All condition params via Michael's email
        %params =[0.62416337 0.001524063 0.69507251]; %prop_spread  tradeoff disrim
        opt_params =[0.06442844 0.1209529 37.93337];
        names={'prop_spread'; 'tradeoff'; 'discrim'};
    case 'qlearning'
        %All condition params via Michael's email
        opt_params =[0.98763,0.216741,0.1448303,0.9854957]; %gamma, alpha, epsilon lambda
        search_params = [.99 fitted_vars.subj_fitting.qlearning.best_parameters(sub,1) fitted_vars.subj_fitting.qlearning.best_parameters(sub,2) .99];
        names={'gamma'; 'alpha'; 'epsilon'; 'lambda'};
    case 'sarsa'
        %All condition params via Michael's email
        %params =[0.989816,0.239363,0.260227,0.9453975]; %gamma, alpha, epsilon lambda
        opt_params =[0.239363,0.260227]; %gamma, alpha, epsilon lambda
        names={'alpha'; 'epsilon';};
    otherwise
        error('Not any agent I''ve heard of');
end
