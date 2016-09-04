function [totalSqErr, ret_all] = TC_minSE(params, blockData, model)

%adaptation of MJF's original code to fit multiple runs of a single subject
%blockData is a struct array containing k elements, where k is number of blocks

totalSqErr=0;

%Identify the number of sessions/blocks

%how many blocks do we have to be fitted
blocks = length(blockData);

ret_all = [];

%for each subject, reset expected value, go, and no-go to zero
priors = [];
priors.V = 0;
priors.Go = 0;
priors.NoGo = 0;

%grand mean RT across all blocks
allRTs = [blockData.rt];
avgRT = mean(allRTs(:));

%fit each block
for b = 1:blocks
    rts = blockData(b).rt; %rt vector
    reward = blockData(b).rew;
    rewFunc = blockData(b).cond; %reward condition
    emo = 'unused';
        
    [RTpred, ret] = TC_Alg(rts, reward, params, priors, avgRT, rewFunc, emo, model);
    
    ret.block = b;
    
    % When fitting multiple blocks within a given subject, use expected value from last trial
    % of block t as the expected value of the first trial for block t + 1.
    priors.V = ret.ev(end);
    
    ret_all = [ret_all; ret];
    
    RTSqErr = (rts - RTpred).^2;
    
    totalSqErr = totalSqErr + sum(RTSqErr);
end