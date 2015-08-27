%rescale the truncated basis so that the typical RBF (not truncated) has range 0..1. The current
%implementation is leading to tiny eligibility traces once you multiply by learning rate.
%technically, this needs to be computed for a given basis setup since the AUC=1.0 scaling will depend on the
%number of timesteps. This correction still maintains greater weighting at the edge (i.e., it is proportionate
%to the AUC=1.0) while giving the reasonable behavior of eligibility representing a "proportion of information
%transfer possible" interpretation.

%find the center that is closest to the midpoint of the time interval
%use the value of the regular basis (0..1 scaling) divided by the truncated basis at the max as the adjustment factor
%note that maxauc above will be equal to this correction factor for all RBFs that are not truncated! Thus, we
%are "undoing" the rescaling of these RBFs, while maintaining the correction at the edge.
%[~,indmid]=min(abs(c - median(tvec)));
%trunc_adjust=max(gaussmat(indmid,:))/max(gaussmat_trunc(indmid,:));
%gaussmat_trunc=gaussmat_trunc*trunc_adjust;


%NB: The working combination is: Regular Gaussian RBF for function evaluation, but weight update proceed by truncated
%Gaussian basis and truncated Gaussian spread function. Conceptually, this is like squeezing in all of the variance of
%the update inside the finite interval but then evaluating the function on the regular RBF... to be developed further
%since this remains a bit of a mystery. :)

%note: with truncated gaussian basis (which functions properly), the idea of extending the basis representation outside
%of the finite interval of interest does not apply (i.e., the functions are not defined outside of the bounds).
%gauss mat for all possible timesteps
%lowest_t = tmin - 4*sig; %3 SDs below lowest center.
%highest_t = tmax + 4*sig;

%t_all = round(lowest_t):round(highest_t);

%gaussmat_all = zeros(nbasis,length(t_all));

%for j = 1:nbasis
%    gaussmat_all(j,:) = gaussmf(t_all,[sig c(j)]);
%end
%plot(t_all, gaussmat_all);


%use this to rescale elibility below to maintain constant AUC equivalent to a standard Gaussian membership function
%this leads to eligibility 0-1.0 for eligibility functions within the interval, and > 1.0 max for truncated functions.
%in testing (weightfit.m), this gives the best sampling behavior


%10/16/2014 NB: Even though some matrices such as v_it have columns for discrete timesteps, the
%agent now learns entirely on a continuous time basis by maximizing the value and uncertainty curves over the
%radial basis set and estimating total uncertainty as the definite integral of the uncertainty function over the
%time window of the trial.
