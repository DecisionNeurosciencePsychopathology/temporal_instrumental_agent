files = glob('/Volumes/bek/vba_results/SHIFTED_U*fixed_uv*.mat');

L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    tau(i) = posterior.muTheta(1)./1000;
end