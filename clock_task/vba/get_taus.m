files = glob('E:/data/sceptic/vba_out/new_lambda_results/SHIFTED_U*kalman_uv*exponential_autocorrel*.mat');

L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    tau(i) = posterior.muTheta(1)./1000;
end


files = glob('E:/data/sceptic/vba_out/new_lambda_results/SHIFTED_U*fixed_decay*exponential_autocorrel*.mat');

L = [];
gamma = [];
for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    chi(i) = posterior.muPhi(end)./100;
end