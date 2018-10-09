function [theta_trans] = transform_theta(theta, inF)

%for now, transforming all elements of theta... should eventually trap models as below, but a bit of a pain
theta_trans = sigmoid(theta); %sigmoid transform alpha and gamma, if they exist

% if any(strcmpi(inF.model, {'fixed', 'decay', 'fixed_UV', 'decay_factorize', 'decay_uniform', ...
% 			     'decay_ps_equate', 'decay_uniform_ps_equate' }))
%   theta_trans = sigmoid(theta); %sigmoid transform alpha, which is theta(1)
% else
%   error(['unrecognized model in transform_theta: ', inF.model]);
% end

end
