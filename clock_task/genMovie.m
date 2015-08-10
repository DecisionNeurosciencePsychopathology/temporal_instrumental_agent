%Wrapper script to make all movies
load('s.mat');
agents = fieldnames(s);
conds = {'IEV' 'DEV' 'QUADUP'};
reversal = 1;
for i = 1:length(agents)
    for j = 1 %1:length(conds)
        clf;
        if reversal
            name = [agents{i} '_' conds{j} '_Reversal'];
        else
            name = [agents{i} '_' conds{j}];
        end
        fprintf('Making movie for %s', name, ' ');
        makeModelMovie(s,agents{i},conds{j},name,reversal);
    end
end

close all;