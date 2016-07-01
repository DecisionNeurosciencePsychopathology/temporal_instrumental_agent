function [c, sig, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(ntimesteps, nbasis, prop_spread, margin_offset, basis_overlap, plot_rbf)
% Shared function for setting up radial basis functions.
%
% Inputs:
%    ntimesteps:    integer representing number of timesteps to be represented in RBF (typically 400)
%    nbasis:        number of basis functions spread evenly over interval (default: 24)
%    prop_spread:   SD of temporal spread function as a proportion of discrete time interval (default: .05)
%    margin_offset: extend the basis interval above and below min and max times by this proportion (default: .125)
%    basis_overlap: Overlap of adjacent RBFs in SD units (Cohen's d).
%    plot_rbf:      Whether to plot RBFs (default: 0)
%
% Outputs:
%    c:             vector of RBF centers in timestep units
%    tvec:          vector of timesteps (1..ntimesteps)
%    sig_spread:    SD of temporal spread function in units of the interval represented (as opposed to proportion)
%    refspread:     Area under the curve of a non-truncated temporal spread function used to rescale AUC of spread to 1.0
%
% Setup centers (means) and sds of basis functions
% Based on testing in fix_rbf_basis.m, default is to place the lowest center 12.5% below the first timestep
% and the last center 12.5% above last timestep. SD should be calculated to give a Cohen's d of
% 1.52 between basis functions (~45% distribution overlap).

if nargin < 2, nbasis=24; end
if nargin < 3, prop_spread=.05; end %default 5% of the interval
if nargin < 4, margin_offset=.125; end %12.5% extension beyond min/max times
if nargin < 5, basis_overlap=1.52; end %default of cohen's d of 1.52 between basis functions
if nargin < 6, plot_rbf=0; end

%Initialize time step vector and allocate for memory
tvec=1:ntimesteps;
sig_spread=prop_spread*range(tvec); %determine SD of spread function

margin_offset = (max(tvec) - min(tvec)).*margin_offset; % convert margin_offset into time scale of tvec

%define lowest and highest centers
tmin = min(tvec) - margin_offset; tmax=max(tvec) + margin_offset;
c=tmin:(tmax-tmin)/(nbasis-1):tmax;

sig = (c(2) - c(1))/basis_overlap;

%construct radial basis matrix using Gaussians
gaussmat = zeros(nbasis, length(tvec));

for j = 1:nbasis
    gaussmat(j,:) = gaussmf(tvec,[sig c(j)]);
end

%version of gaussian where each function has AUC = 1.0 (PDF representation). Not currently used (MH Sep2015)
maxauc_all=max(sum(gaussmat, 2));
%gaussmat_pdf=gaussmat./maxauc_all;

%normalize RBFs to each have AUC = 1.0 within observed time interval
%this is essentially a truncated Gaussian basis such that AUC = 1.0 for all basis functions within the discrete interval
maxauc_each=sum(gaussmat,2)*ones(1,length(tvec)); %outer product of vectors to allow for col-wise division below
gaussmat_trunc=gaussmat./maxauc_each;
%gaussmat_trunc=gaussmat_pdf;

%basis functions at the edge are too sharp -- the agent is inflating value toward infinity. Need to soften the tails.
%what if it's AUC=1.0 divided by some proportion of basis in the interval such that less gets downweighted?

if plot_rbf == 1
    figure(20); plot(tvec,gaussmat); title('Regular RBF');
    figure(21); plot(tvec,gaussmat_trunc); title('Truncated RBF');
end

%Determine the AUC of a non-truncated eligilibity function
%Use this to rescale elibility function to maintain constant AUC equivalent to a standard Gaussian membership function.
%This leads to eligibility 0-1.0 for eligibility functions within the interval, and > 1.0 max for truncated functions.
%In testing (weightfit.m), this gives the best sampling behavior
refspread = sum(gaussmf(min(tvec)-range(tvec):max(tvec)+range(tvec), [sig_spread, median(tvec)]));

end