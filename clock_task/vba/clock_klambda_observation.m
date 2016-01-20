function  [ gx ] = clock_tc_observation(hidden,phi,u,inG)
% K+Lambda only observation (choice) function

% INPUT
%   - hidden : evolved hidden states from evolution function: NONE
%
%   - phi : observation parameters
%       (1) K (intercept)
%       (2) lambda (autocorrelation)
%       
%   - u (inputs) :
%       (1) RT (RT on trial t)
%       (2) RT_prev (RT[t-1] -- used in observation function for resetting explore)
%
% - inG : input structure
%
% OUTPUT
% - gx : predicted RT

%pull out input parameters
%RT = u(1); %current RT
RT_prev = u(2); %previous RT

K = unifinv(fastnormcdf(phi(1)), 0, inG.maxRT);
lambda = 1 ./ (1+exp(-phi(2))); %exponential transform to 0..1

gx = K + lambda*RT_prev;