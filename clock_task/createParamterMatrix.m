function data=createParamterMatrix
%Script to grab every variable from vba output and put it into a struct

 modelnames = {'qlearning','fixed','fixed_decay', 'kalman_softmax', 'kalman_processnoise'...
    'kalman_logistic','kalman_sigmavolatility','kalman_uv_sum','fixed_uv'};

%  modelnames = {'fixed','fixed_decay', 'kalman_softmax', 'kalman_processnoise'...
%     'kalman_logistic','kalman_sigmavolatility','kalman_uv_sum','fixed_uv'};
% % modelnames = {'kalman_logistic','kalman_sigmavolatility','kalman_uv_sum','fixed_uv','qlearning'};
%modelnames = {'kalman_processnoise'}; %Forgot omega

%Grab everyone but qlearning and kalman_processnoise for now
%modelnames = {'fixed','fixed_decay', 'kalman_softmax','kalman_processnoise',...
 %   'kalman_logistic','kalman_sigmavolatility','kalman_uv_sum','fixed_uv',};

%arg_string = '/SHIFTED_U_CORRECT*_multinomial1_multi*0_fixed*1_ua*1_sc*fixed_prop_spread.mat';

data = struct;

for i = 1:length(modelnames)
    
    arg_string = sprintf('/SHIFTED_U_CORRECT*_%s_multinomial1_multi*0_fixed*1_ua*1_sc*fixed_prop_spread.mat',modelnames{i}); %Not organzied by model

    
    %Initialize storage vars
    evo_params = [];
    obs_params = [];
    
    %Grab the files
    %fnames = glob(['/Volumes/bek/vba_results/', modelnames{i},arg_string]);
    fnames = glob(['/Volumes/bek/vba_results/', arg_string]); %Not organzied by model
    
    %Special qlearning case was called q_step
    if strcmp(modelnames{i},'qlearning')
        arg_string = '/SHIFTED_U_CORRECT*_multinomial1_multi*1_fixed*1*.mat';
        fnames = glob(['/Volumes/bek/vba_results/',modelnames{i}, arg_string]); %Not organzied by model
    end
    

    
    %Load the files
    fprintf('\nLoading the files and gathering the parameters')
    for j = 1:length(fnames)
        load(fnames{j})
        
        %Cbeck if propspread was used, it should have been, if not break
        %for now...
        %06/14/2016 we decided to just fix propspread moving forward
%         if out.options.inF.fit_propspread~=1
%             return
%         end
        
        %Grab the needed data
        evo_params = [evo_params posterior.muTheta];
        obs_params = [obs_params posterior.muPhi];
        fprintf('.')
    end
    fprintf('Done!\n\n')
    data = organizeData(data,modelnames{i},evo_params,obs_params,out);
    
end


function s=organizeData(s,name,theta,P,out)
%Set the proper paramters NOTE gamma is first then phi for sigma volatility
%model.

%We could probably just add the beta for everyone up here because it
%appears it is P(1) for everyone...but just to be sure.


%If we are using qlearning
if strcmp(name,'qlearning')
    s.transformed.(name).alpha = 1./(1+exp(-theta(1)));
    s.non_transformed.(name).alpha = theta(1);
    s.transformed.(name).beta = exp(P);
    s.non_transformed.(name).beta = P;
else
    
    %KalmanLogistic / all obs function params
    if out.options.inF.kalman.kalman_logistic
        s.transformed.(name).beta = exp(P(1,:));
        s.non_transformed.(name).beta = P(1,:);
        
        s.transformed.(name).discrim = 1./(1+exp(-P(2,:)));
        s.non_transformed.(name).discrim = P(2,:);
    else
        s.transformed.(name).beta = exp(P);
        s.non_transformed.(name).beta = P;
    end
    %ProcessNoise
    %if out.options.inF.kalman.kalman_processnoise %Correct
    if out.options.inF.kalman.processnoise
        s.transformed.(name).omega = 1./(1+exp(-theta(1,:)));
        s.non_transformed.(name).omega = theta(1,:);
    end
    %SigVol
    if out.options.inF.kalman.kalman_sigmavolatility || out.options.inF.kalman.kalman_sigmavolatility_local || out.options.inF.kalman.kalman_sigmavolatility_precision
        if out.options.inF.no_gamma
            s.(name).gamma = 1-phi;
        else
            s.transformed.(name).gamma = 1./(1+exp(-theta(2,:)));
            s.non_transformed.(name).gamma = theta(2,:);
        end
        s.transformed.(name).phi = 1./(1+exp(-theta(1,:)));
        s.non_transformed.(name).phi = theta(1,:);
    end
    %UVSumSigVol
    if out.options.inF.kalman.kalman_uv_sum_sig_vol
        s.transformed.(name).gamma = 1./(1+exp(-theta(2,:)));
        s.non_transformed.(name).gamma = theta(2,:);
        
        s.transformed.(name).phi = 1./(1+exp(-theta(3,:)));
        s.non_transformed.(name).phi = theta(3,:);
    end
    %UVLogistic
    if out.options.inF.kalman.kalman_uv_logistic
        s.transformed.(name).tradeoff = 1./(1+exp(-theta(1,:)));
        s.non_transformed.(name).tradeoff = theta(1,:);
    end
    %Fixed_uv
    if out.options.inF.kalman.fixed_uv
        s.transformed.(name).alpha = 1./(1+exp(-theta(2,:)));
        s.non_transformed.(name).alpha = theta(2,:);
    end %Alpha should always be two as tau should always be 1
    
    %UVSum
    if out.options.inF.kalman.kalman_uv_sum || out.options.inF.kalman.kalman_uv_sum_sig_vol || out.options.inF.kalman.fixed_uv
        if out.options.inF.u_aversion
            s.transformed.(name).tau = theta(1,:)./1000; % scale it down a bit
            s.non_transformed.(name).tau = theta(1,:); % scale it down a bit
        end
    end
    %FixedDecay
    if strcmp(name,'fixed_decay')
        s.transformed.(name).alpha = 1./(1+exp(-theta(1,:)));
        s.non_transformed.(name).alpha = theta(1,:);
        
        s.transformed.(name).gamma = (1./(1+exp(-theta(2,:))))./1;
        s.non_transformed.(name).gamma = theta(2,:);
    end
    %Fixed
    if strcmp(name,'fixed')
        s.transformed.(name).alpha = 1./(1+exp(-theta(1,:)));
        s.non_transformed.(name).alpha = theta(1,:);
    end
    %All propspread
    last_idx = size(theta,1);
    if out.options.inF.fit_propspread
        s.transformed.(name).prop_spread = 1./(1+exp(-theta(last_idx,:)));
        s.non_transformed.(name).prop_spread = theta(last_idx,:);
    else
        s.(name).prop_spread = out.options.inF.sig_spread;
    end
end