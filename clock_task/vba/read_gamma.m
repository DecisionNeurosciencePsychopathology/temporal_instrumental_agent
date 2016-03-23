files = glob('/Volumes/bek/vba_results/SHIFTED_U*fixed_decay*.mat');

L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    gamma(i) = (1./(1+exp(-posterior.muTheta(2)))); 
    
end