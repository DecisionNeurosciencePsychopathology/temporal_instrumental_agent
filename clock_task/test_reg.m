
function s = test_reg(rts_obs,ret,modelname, step)
%% takes in actual rts, predicted RTexploit and predicted RTexplore, 
% fits regression, prints out results 
% returns struct s with
% s.b s.bint s.r s.rint s.stats where
% s.b = betas
% s.bint = beta CIs
% s.r = residuals
% s.rint = residual CIs
% s.stats = [R2 F-stat p-value estimate_of_the_error_variance]
% can do stepwise regression if step=1
if nargin<4
    step = 0;
end



%%
if strcmpi(modelname,'value_softmax')
X = [ret.rts_pred_exploit(1:50)' ones(size(ret.rts_pred_explore(1:50)'))];
else
X = [ret.rts_pred_exploit(1:50)' ret.rts_pred_explore(1:50)' ones(size(ret.rts_pred_explore(1:50)'))];
end

y = round(rts_obs./10);
[s.b, s.bint, s.r, s.rint, s.stats] = regress(y, X);

if strcmpi(modelname,'value_softmax')
fprintf('exploit_beta = %.2f,  CI=%.2f - %.2f \r',s.b(1), s.bint(1,1), s.bint(1,2));
fprintf('R2 = %.3f,  F=%.2f,  p = %.5f\r\r',s.stats(1), s.stats(2), s.stats(3))
else
fprintf('exploit_beta = %.2f,  CI=%.2f - %.2f; explore_beta = %.2f CI=%.2f - %.2f\r',s.b(1), s.bint(1,1), s.bint(1,2), s.b(2), s.bint(2,1), s.bint(2,2));
fprintf('R2 = %.3f,  F=%.2f,  p = %.5f\r\r',s.stats(1), s.stats(2), s.stats(3))
end

%%
if step == 1
stepwise(X,y);
end
end