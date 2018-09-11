function [theta_trans] = transform_theta(theta, inF)

if any(strcmpi(inF.model, {'fixed', 'decay', 'fixed_UV', 'decay_factorize'}))
  theta_trans = sigmoid(theta); %sigmoid transform alpha, which is theta(1)
else
  error(['unrecognized model in transform_theta: ', inF.model]);
end

%handle other models at some point...

end
