function [phi_trans] = transform_phi(phi, inG)

%force positive temperature (that's the only element of phi)
phi_trans = exp(phi);

%handle other possibilities at some point...

end