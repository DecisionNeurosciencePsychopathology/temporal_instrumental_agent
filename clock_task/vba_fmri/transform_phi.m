function [phi_trans] = transform_phi(phi, inG)

if ismember(inG.model, {'fixed_uv', 'fixed_uv_ureset', 'fixed_uv_ureset_fixedparams_fmri', 'fixed_uv_ureset_fixedparams_meg'})
  phi_trans = [ exp(phi(1)); ... %temperature
    phi(2)./100 ]; % scale down tau since it tends to be << 1, but in fitting, having a larger range may help estimation
else
  %force positive temperature (that's the only element of phi)
  phi_trans = exp(phi);
end


%handle other possibilities at some point...

end
