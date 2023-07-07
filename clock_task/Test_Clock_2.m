% 2023-04-12 AndyP
% Test N Mushrooms Clock 2.0 vs N Points Clock 1.0

condset = {'Inquisit_IEV','Inquisit_DEV','Inquisit_CEV','Inquisit_CEVR'};
nSims = 500;
RT = linspace(0,5000,50);
%rt = repmat(RT',[length(condset),1]);
%rt = repmat(rt,[nSims,1]);

rew = [];
ev = [];
frq = [];
mag = [];
sim = [];
mush = [];
cond = [];
rt = [];

doInquisit_POS = 1;
gamma = 5000/100;

for iS = 1:nSims
    disp(iS)
    for iC=1:length(condset)
        cond0 = condset{iC};
        %disp(cond0);
        rng('shuffle');
        for iT=1:length(RT)
            
            if doInquisit_POS==1
                if mod(2,iS)==0
                    start_pos = 0 + (75-0)*rand(1,1);
                else
                    start_pos = 75 + (100-75)*rand(1,1);
                end
                if start_pos >=0 && start_pos < 75
                    RT(iT) = mod(gamma*(start_pos+25)+RT(iT),5000);
                else
                    RT(iT) = mod(gamma*(start_pos-75)+RT(iT),5000);
                end
            end
            
            [rew0, ev0, frq0, mag0] = RewFunction1(RT(iT), cond0, 0, 5000);
            rew = cat(1,rew,rew0);
            ev = cat(1,ev,ev0);
            frq = cat(1,frq,frq0);
            mag = cat(1,mag,mag0);
            sim = cat(1,sim,iS);
            rt = cat(1,rt,RT(iT));
            
            if strcmp(cond0,'IEV') || strcmp(cond0,'Inquisit_IEV') || strcmp(cond0,'Explore_IEV')
                cond1 = 1;
            elseif strcmp(cond0, 'DEV') || strcmp(cond0,'Inquisit_DEV') || strcmp(cond0,'Explore_DEV')
                cond1 = 2;
            elseif strcmp(cond0, 'CEV') || strcmp(cond0,'Inquisit_CEV') || strcmp(cond0,'Explore_CEV')
                cond1 = 3;
            elseif strcmp(cond0, 'CEVR') || strcmp(cond0,'Inquisit_CEVR') || strcmp(cond0,'Explore_CEVR')
                cond1 = 4;
            end
            cond = cat(1,cond,cond1);
            
            
            if rew0 <= 24
                mush0 = 1;
            elseif rew0 > 24 && rew0 <= 32
                mush0 = 2;
            elseif rew0 > 32 && rew0 <= 46
                mush0 = 3;
            elseif rew0 > 46 && rew0 <=60
                mush0 = 4;
            elseif rew0 > 60 && rew0 <= 74
                mush0 = 5;
            elseif rew0 > 74 && rew0 <= 88
                mush0 = 6;
            elseif rew0 > 88 && rew0 <= 102
                mush0 = 7;
            elseif rew0 > 102 && rew0 <= 124
                mush0 = 8;
            elseif rew0 > 124 && rew0 <= 136
                mush0 = 9;
            elseif rew0 > 136
                mush0 = 10;
            end
            
            mush = cat(1,mush,mush0);
        end
    end
end