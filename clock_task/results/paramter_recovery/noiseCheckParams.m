std_mean_error=std(fitted_vars.noise_params.fixed_rho.best_parameters)/sqrt(100);

model = 'v_processnoise';
%org_params = [.2 .002 .0007];
org_params = 10;
nosie_sse=sqrt(sum((repmat(org_params,100,1)-fitted_vars.noise_params.(model).best_parameters).^2)./100);
%sqrt(sum((0.2-fitted_vars.noise_params.(model).best_parameters(:,1)).^2)/100)

%Look at params via hist
figure(88)
hist(fitted_vars.noise_params.(model).best_parameters(:,1));
figure(89)
hist(fitted_vars.noise_params.(model).best_parameters(:,2));
figure(90)
hist(fitted_vars.noise_params.(model).best_parameters(:,3));