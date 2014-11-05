function cost = make_gaussmat(params)
%% Makes a gaussian basis function adding up to zero or 1
% the four params are the amplitudes of the two first and two last gaussian
ntimesteps = 2000;
t=1:ntimesteps;
sig = 150;
number_of_stimuli = 12;
step_num = number_of_stimuli - 1;
step = (ntimesteps-1)/step_num;
c = -step:step:ntimesteps+step;
nbasis=length(c);
gaussmat = zeros(nbasis,ntimesteps);

for j = 1:nbasis
    gaussmat(j,:) = gaussmf(t,[sig c(j)]);
end

    temp_vh = ones(1,ntimesteps);
    temp_vh = [params(1:2) ones(1,10) params(3:4)]'*temp_vh;
corrected_basis=temp_vh.*gaussmat; 
resulting_curve = sum(corrected_basis);

%% the idea is to bring the edges up to the max of sum(gaussmat)
curve_discrepancy = abs(resulting_curve-max(sum(gaussmat)));
cost = sum(curve_discrepancy);