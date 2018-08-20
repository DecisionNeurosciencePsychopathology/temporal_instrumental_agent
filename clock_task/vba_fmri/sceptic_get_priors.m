function [priors] = sceptic_get_priors(dim)
  priors=[];
  priors.a_alpha = Inf;   % infinite precision prior
  priors.b_alpha = 0;
  priors.a_sigma = 1;     % Jeffrey's prior
  priors.b_sigma = 1;     % Jeffrey's prior
  
  priors.muPhi = zeros(dim.n_phi,1); % exp tranform on temperature inside observation fx
  priors.SigmaPhi = 1e1*eye(dim.n_phi); %variance of 10

  priors.muTheta = zeros(dim.n_theta,1); %mean of 0
  priors.SigmaTheta = 1e1*eye(dim.n_theta); %variance of 10

  %0 mean and variance on initial states
  priors.muX0 = zeros(dim.n,1);
  priors.SigmaX0 = zeros(dim.n);

end
