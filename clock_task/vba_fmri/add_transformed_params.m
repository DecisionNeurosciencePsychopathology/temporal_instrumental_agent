function [posterior] = add_transformed_params(posterior, so)

posterior.transformed.muPhi = transform_phi(posterior.muPhi, so);
posterior.transformed.SigmaPhi = transform_phi(posterior.SigmaPhi, so);

posterior.transformed.muTheta = transform_theta(posterior.muTheta, so);
posterior.transformed.SigmaTheta = transform_theta(posterior.SigmaTheta, so);

%in MFX case, populate transformed ffx variants
if isfield(posterior, 'ffx')
  posterior.transformed.muPhi_ffx = transform_phi(posterior.ffx.muPhi, so);
  posterior.transformed.SigmaPhi_ffx = transform_phi(posterior.ffx.SigmaPhi, so);
  
  posterior.transformed.muTheta_ffx = transform_theta(posterior.ffx.muTheta, so);
  posterior.transformed.SigmaTheta_ffx = transform_theta(posterior.ffx.SigmaTheta, so);
end

end