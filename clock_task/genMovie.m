%Wrapper script to make all movies
load('s.mat');
agents = fieldnames(s);
conds = {'IEV' 'DEV' 'QUADUP'};
for i = 4:length(agents)
    for j = 1:length(conds)
        clf;
        name = [agents{i} '_' conds{j}];
        fprintf('Making movie for %s', name);
        makeModelMovie(s,agents{i},conds{j},name)
    end
end