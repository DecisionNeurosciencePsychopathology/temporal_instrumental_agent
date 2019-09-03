function [theta_trans] = transform_theta(theta, inF)

%only fixed uv baked has tau in theta (part of learning rule)
if ismember(inF.model, {'fixed_uv_baked', 'fixed_uv_baked_ureset'})
  %fixed_uv is the variant where tau is allowed to be positive or negative
  %fixed_uv_positive would be the version where we enforce positiviety

  theta_trans = [ VBA_sigmoid(theta(1)); ... %alpha
    theta(2)./100 ]; % scale down tau since it tends to be << 1, but in fitting, having a larger range may help estimation
elseif strcmpi(inF.model, 'fixed_uv_positive')
  %this is just a stub in the port; I have not implemented this in full
  theta_trans = [ VBA_sigmoid(theta(1)); ... %alpha
    1./(1+exp(-theta(2)-log(inF.sigma_noise))) ]; %tau -- discounted by noise?
else
  %for now, transforming all elements of theta... should eventually trap models as below, but a bit of a pain
  theta_trans = VBA_sigmoid(theta); %sigmoid transform alpha and gamma, if they exist
end



% if any(strcmpi(inF.model, {'fixed', 'decay', 'fixed_UV', 'decay_factorize', 'decay_uniform', ...
% 			     'decay_ps_equate', 'decay_uniform_ps_equate' }))
%   theta_trans = sigmoid(theta); %sigmoid transform alpha, which is theta(1)
% else
%   error(['unrecognized model in transform_theta: ', inF.model]);
% end

end
