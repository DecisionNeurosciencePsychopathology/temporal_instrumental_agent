std_mean_error=std(fitted_vars.noise_params.fixed_rho.best_parameters)/sqrt(100);

agent = 'fixedLR_softmax'; %change angent name here
%org_params = [.2 .002 .0007];
org_params = [0.03081894 0.01698200];
len = length(fitted_vars.noise_params.(agent).best_parameters);
nosie_sse=sqrt(sum((repmat(org_params,len,1)-fitted_vars.noise_params.(agent).best_parameters).^2)./100);
%sqrt(sum((0.2-fitted_vars.noise_params.(agent).best_parameters(:,1)).^2)/100)

%Look at params via hist
figure(88)
hist(fitted_vars.noise_params.(agent).best_parameters(:,1));
figure(89)
hist(fitted_vars.noise_params.(agent).best_parameters(:,2));
figure(90)
hist(fitted_vars.noise_params.(agent).best_parameters(:,3));


for j = 1:size(fitted_vars.noise_params.(agent).best_parameters,2)
    figure(10+j)
    hist(fitted_vars.noise_params.(agent).best_parameters(:,j),20);
end