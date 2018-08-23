function [theta_trans] = transform_theta(theta, inF)

if strcmpi(inF.model, 'fixed')
  theta_trans = sigmoid(theta); %sigmoid transform alpha, which is theta(1)
elseif strcmpi(inF.model, 'decay')
  theta_trans = sigmoid(theta); %sigmoid transform alpha and gamma
end

%handle other models at some point...

end