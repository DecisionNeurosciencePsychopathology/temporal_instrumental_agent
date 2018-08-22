function [posterior] = add_transformed_params(posterior, so)

posterior.transformed.muPhi = transform_phi(posterior.muPhi, so);
posterior.transformed.SigmaPhi = transform_phi(posterior.SigmaPhi, so);

posterior.transformed.muTheta = transform_theta(posterior.muTheta, so);
posterior.transformed.SigmaTheta = transform_theta(posterior.SigmaTheta, so);

end