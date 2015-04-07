%MNH Mar2015: figure out how to approximate uniform random sampling by moving to point of greatest uncertainty on each
%trial. The goal is for the sampling strategy to approximate a random uniform distribution over time, but with the
%sequential sampling always moving to point of greatest uncertainty. The major problem so far is that even with a good
%radial basis set (see fix_rbf_basis.m for how the basis was perfected to fit any number of functions), the process of
%propagating sampling/value information throughout the finite time interval leads to highly uneven sampling with a
%tendency to oversample the edges due to insufficient weight updates at the edge.

%here, I demonstrate the basic problem of an even weighting scheme in a finite interval RBF basis. This echos papers
%from John Boyd demonstrating problems with RBF fitting over a finite interval. One initial approach (below) was to look
%at the RBF weighting scheme needed to approximate f(x) === 1.0 and to use these values as a correction factor for
%weight updates. Although this is a promising approach when trying to fit a known function (e.g., f(x) === 1.0), it does
%not perform well with sequential sampling (i.e., with a trial-wise temporal spread function).

%instead, the successful strategy is to define a radial basis where each function has equal area under the curve within
%the observed time interval (here, AUC=1.0 for all RBFs). Conceptually, this means that an equal weighting scheme should
%generate a relatively flat function approximation even at the edge because an equivalent weight acting over a truncated
%function will lead to an effectively higher height/contribution in the RBF summation. Likewise, to propagate value
%information evenly over the finite interval, a truncated spread function always having AUC=1.0 is necessary to evenly
%distribute the sampling information within the interval.

%there are nuances to be worked out in terms of number of basis functions, degree of added Gaussian noise, separation
%between basis functions, etc., but the truncated Gaussian RBF with truncated Gaussian spread function works well for
%many parameter sets (esp. as long as the number of basis function is not low, such as 8).


%%
%demonstration of basic problem with rbf with equal weights...
%even with large offsets, there will always be falloff at the edges
crazycenters=-150:50:650;
res=rbfeval(1:500, ones(1,length(crazycenters)), crazycenters, ones(1,length(crazycenters)).*66);
plot(res)

%falloff persists until the entire symmetric interval about the min and max is covered (i.e., 100% extension
%on each side). See below: this looks fine
crazycenters=-500:50:1000;
res=rbfeval(1:500, ones(1,length(crazycenters)), crazycenters, ones(1,length(crazycenters)).*66);
plot(res)


%try with truncated gaussians
%get an alternate weird pattern with greater heights at the edges
crazycenters=0:50:500;
res=rbfeval(1:500, ones(1,length(crazycenters)), crazycenters, ones(1,length(crazycenters)).*33, [0 500]);
plot(res)

%check the value of having at least one basis function outside the observed interval
spacing=38.461535;
sig=25.64; %cohen's d = 1.5
centers_0offset=0:spacing:500;
centers_1offset=-spacing:spacing:500+spacing; %14 even basis functions with no offset
centers_2offset=-spacing*2:spacing:500+spacing*2; %14 even basis functions with no offset
%centers=0:spacing:500; %14 even basis functions with no offset
%centers=0:25:500;
init_params_0=ones(1,length(centers_0offset));
upper_bounds_0=ones(1,length(centers_0offset)).*10;
lower_bounds_0=ones(1,length(centers_0offset)).*0;

init_params_1=ones(1,length(centers_1offset));
upper_bounds_1=ones(1,length(centers_1offset)).*10;
lower_bounds_1=ones(1,length(centers_1offset)).*0;

init_params_2=ones(1,length(centers_2offset));
upper_bounds_2=ones(1,length(centers_2offset)).*10;
lower_bounds_2=ones(1,length(centers_2offset)).*0;

fmincon_options = optimoptions(@fmincon,'TolX',1e-10,'TolFun',1e-10,'MaxIter',100000, 'MaxFunEvals', 100000); %request a bit higher precision (doesn't really matter here!)
%fmincon_default = optimoptions(@fmincon);

%[par1, cost1, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers, ones(1,length(centers)).*sig, -Inf), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_default);
[par_0offset, cost_0offset, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_0offset, ones(1,length(centers_0offset)).*sig, -Inf), init_params_0, [], [], [], [], lower_bounds_0, upper_bounds_0, [], fmincon_options);
[par_1offset, cost_1offset, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_1offset, ones(1,length(centers_1offset)).*sig, -Inf), init_params_1, [], [], [], [], lower_bounds_1, upper_bounds_1, [], fmincon_options);
[par_2offset, cost_2offset, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_2offset, ones(1,length(centers_2offset)).*sig, -Inf), init_params_2, [], [], [], [], lower_bounds_2, upper_bounds_2, [], fmincon_options);

%cost for 1 basis function outside of observed interval is noticeably better, but 2 is only a trivial
%improvement

f0 = rbfeval(1:500, par_0offset, centers_0offset, sig, -Inf);
f1 = rbfeval(1:500, par_1offset, centers_1offset, sig, -Inf);
f2 = rbfeval(1:500, par_2offset, centers_2offset, sig, -Inf);
figure(1); plot(f0, 'b'); hold on
plot(f1, 'g'); 
plot(f2, 'r'); hold off;

%check corrective basis using truncated gaussians
init_params_0=ones(1,length(centers_0offset));
upper_bounds_0=ones(1,length(centers_0offset)).*50;
lower_bounds_0=ones(1,length(centers_0offset)).*0;

init_params_1=ones(1,length(centers_1offset));
upper_bounds_1=ones(1,length(centers_1offset)).*50;
lower_bounds_1=ones(1,length(centers_1offset)).*0;

init_params_2=ones(1,length(centers_2offset));
upper_bounds_2=ones(1,length(centers_2offset)).*50;
lower_bounds_2=ones(1,length(centers_2offset)).*0;

[par_0offset_trunc, cost_0offset_trunc, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_0offset, ones(1,length(centers_0offset)).*sig, [1 500]), init_params_0, [], [], [], [], lower_bounds_0, upper_bounds_0, [], fmincon_options);
[par_1offset_trunc, cost_1offset_trunc, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_1offset, ones(1,length(centers_1offset)).*sig, [1 500]), init_params_1, [], [], [], [], lower_bounds_1, upper_bounds_1, [], fmincon_options);
[par_2offset_trunc, cost_2offset_trunc, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_2offset, ones(1,length(centers_2offset)).*sig, [1 500]), init_params_2, [], [], [], [], lower_bounds_2, upper_bounds_2, [], fmincon_options);

par_0offset_trunc;
par_1offset_trunc
par_2offset_trunc

cost_0offset_trunc
cost_1offset_trunc
cost_2offset_trunc

f0 = rbfeval(1:500, par_0offset_trunc, centers_0offset, ones(1,length(centers_0offset)).*sig, [1 500]);
f1 = rbfeval(1:500, par_1offset_trunc, centers_1offset, ones(1,length(centers_1offset)).*sig, [1 500]);
f2 = rbfeval(1:500, par_2offset_trunc, centers_2offset, ones(1,length(centers_2offset)).*sig, [1 500]);
figure(2); plot(f0, 'b'); hold on
plot(f1, 'g'); 
plot(f2, 'r'); hold off;


%conclusions: fitting a correction to the weights allows for an f(x)===1.0 with reasonable precision over a
%finite interval. Using a truncated Gaussian basis leads to a different set of weights, but the same qualitative
%pattern of fit (and trivial cost differences).
%%

%this is the winner... truncated gaussian basis with truncated gaussian temporal spread function
%so, proportion of temporal receptive field for a given RBF intersected with AUC of the spread function
%and random noise does help -- makes sense in overcoming tiny troughs in the basis function set that will be oversampled
%by deterministic sampling

nbasis=14;
ntrials=1000;
tvec=1:500;

[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 30], nbasis, 10000, 'gauss_auc', tvec, 0, 0, 30);
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 30], nbasis, 1000, 'gauss_auc', tvec, 0, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 30], nbasis, 1000, 'gauss_auc', tvec, 0, 0, 30);

[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .06 .04], nbasis, 500, 'gauss_auc', tvec, 0, 0);

%does the shape of rt sampling without noise depend on nbasis etc?
%a major problem is if the basis functions are too sparse (i.e., not enough), then we cannot recover
%a reasonably flat sampling curve -- the updates are too lumpy
[cost, rts, w, w_update, u_all] = weightfit([0 2.0 .02 0], 8, 100, 'gauss_auc', tvec, 0, 1);
figure(1); hist(rts)
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .06 0], 14, 500, 'gauss_auc', tvec, 0, 0);
figure(2); hist(rts)
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .06 0], 30, 500, 'gauss_auc', tvec, 0, 0);
figure(3); hist(rts)

%having a relatively precise temporal generalization function, even without noise, does well
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .01 0], 30, 5000, 'gauss_auc', tvec, 0, 0);
figure(3); hist(rts)


%having a relatively precise temporal generalization function, even without noise, does well if there are enough basis functions
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .01 0], 30, 500, 'gauss_mult', tvec, 0, 0); %cost=.91
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .01 .02], 30, 500, 'gauss_mult', tvec, 0, 0); %cost=.82


%conclusions: even though I thought the gauss_auc approach was conceptually better, it leads to lower weights
%at the edges and consequent oversampling of the edges... seemingly even with a truncated basis. Need to
%figure this out.




[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .04 .00], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .04 .02], 30, 5000, 'gauss_mult', tvec, 0, 0);
figure(5); hist(rts);

%at these settings, RBFEVAL on AUC1 outperforms RBFEVAL on TRUNC with any level of noise
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .01 .00], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .01 .01], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .01 .02], 30, 5000, 'gauss_mult', tvec, 0, 0);

%at these settings, RBFEVAL on AUC1 outperforms RBFEVAL on TRUNC with any level of noise
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .02 .00], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .02 .005], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .02 .01], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .02 .02], 30, 5000, 'gauss_mult', tvec, 0, 0);

%at these settings, RBFEVAL on AUC1 outperforms RBFEVAL on TRUNC with any level of noise
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .04 .00], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .04 .01], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .04 .02], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .04 .04], 30, 5000, 'gauss_mult', tvec, 0, 0);

%at these settings, RBFEVAL on AUC1 outperforms RBFEVAL on TRUNC with any level of noise
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .08 .00], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .08 .01], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .08 .02], 30, 5000, 'gauss_mult', tvec, 0, 0);
[cost, rts, w, w_update, u_all] = weightfit([.05 1.5 .08 .04], 30, 5000, 'gauss_mult', tvec, 0, 0);




%Apr2015: simulation matrix of potentially good approaches:
%1) gauss_mult or gauss_auc
%2) truncate eligibility to have AUC=1.0 (N.B. actually need to have AUC of a non-truncated RBF Gaussian membership for
%the weight updates to run 0 ~ 1
%3) multiply eligibility by regular or truncated RBF
%4) evaluate weights on regular or truncated RBF

%Incomplete testing above revealed the the winning combination is gauss_mult, truncated eligibility multiplied by
%truncated RBF, then evaluate weights on regular RBF. (So, during weight update, squeeze everything inside observed
%interval). And yes, Gaussian noise is helpful in overcoming uneven sampling! (Often with sigma of 50-100% of the
%temporal generalization function sd).




% GAUSS AUC MAX WEIGHT UPDATES
% ans =
% 
%   Columns 1 through 8
% 
%     0.3794    0.4768    0.6400    0.9181    0.8713    0.8481    0.8439    0.8438
% 
%   Columns 9 through 16
% 
%     0.8438    0.8438    0.8438    0.8438    0.8437    0.8438    0.8438    0.8438
% 
%   Columns 17 through 24
% 
%     0.8438    0.8438    0.8439    0.8481    0.8796    0.9181    0.6400    0.4768
% 
%   Column 25
% 
%     0.3794


% GAUSS MULT MAX WEIGHT UPDATES
% 
% ans =
% 
%   Columns 1 through 8
% 
%     1.9290    1.8921    1.7945    1.5097    0.9383    0.8172    0.8108    0.8106
% 
%   Columns 9 through 16
% 
%     0.8106    0.8106    0.8106    0.8106    0.8106    0.8106    0.8106    0.8106
% 
%   Columns 17 through 24
% 
%     0.8106    0.8106    0.8108    0.8172    0.9383    1.5097    1.7945    1.8921
% 
%   Column 25
% 
%     1.9290


[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .05 .02], 25, 5000, 'gauss_auc', tvec, 3, 0);


%but with a low number of basis functions, absence of noise is a problem, but a little noise helps a lot
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .01 0], 10, 500, 'gauss_auc', tvec, 0, 0); %cost=91
[cost, rts, w, w_update, u_all] = weightfit([.15 1.5 .01 .02], 10, 500, 'gauss_auc', tvec, 0, 0); %cost=26
[cost, rts, w, w_update, u_all] = weightfit([.15 1.4 .01 .02], 10, 500, 'gauss_auc', tvec, 0, 0); %cost=12
[cost, rts, w, w_update, u_all] = weightfit([.15 1.3 .01 .02], 10, 500, 'gauss_auc', tvec, 0, 0); %cost=9.7
[cost, rts, w, w_update, u_all] = weightfit([.15 1.0 .01 .02], 10, 500, 'gauss_auc', tvec, 0, 0); %cost=208 --increase

%some evidence that with a lower number of basis functions, decreasing separation among them (i.e., higher SD)
%is helpful up to a point

%some conclusions: 
%1) noise is helpful, but especially with a low number of basis functions.
%2) having too few basis (<10?) functions is inherently problematic for sampling evenly and representing complex value
%functions
%3) 5-15% offset for the lowest/highest basis functions is sufficient to reduce cost at the edge
%4) separation of about d=1.3 - 1.5 is often useful for a reasonable number of basis functions (15-30)
%5) there is not a huge benefit to adding a lot of basis functions (50, 100, etc.) compared to 15-25
%6) A truncated Gaussian RBF set updated by a truncated Gaussian spread function (where all Gaussians have AUC=1.0)
%yields very reasonable behavior at the edge and throughout the interval.





hist(rts)
%one conclusion to be drawn from this is that noise is crucial to overcoming lumpiness in the basis function if there
%are a rel


%look at optimization
fmincon_options = optimoptions(@fmincon, 'DiffMinChange', 0.01);
fmincon_options.MaxFunEvals = 5000;
fmincon_options.MaxIter = 5000;

lower_bounds=[0, 0.3, .005 .005]; %prop_offset, basis_sep, prop_spread, prop_noise
upper_bounds=[0.5, 3.0 1.0 1.0];
init_params=[0.125, 1.5 .06 .04];


[par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis, ntrials, 'gauss_auc', tvec), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

%params: 0.066998      2.1985   0.0074313    0.037972
%fairly low overlap between functions (leads to lumpy function)
%very precise gaussian update (low spread) also gives funny results...
%what does the agent look like at these params
[cost, rts, w, w_update, u_all]=weightfit(par, nbasis, 120, 'gauss_auc', tvec, 0, 1);

%what about with a more reasonable basis_sep, spread, and noise?
%this agent seems good in terms of temporal precision
[cost, rts, w, w_update, u_all]=weightfit([.07 1.5 .04 .025], nbasis, 120, 'gauss_auc', tvec, 0, 1);




%current promising basis set
%without noise
[cost, rts, w, w_update, u_all]=weightfit([0.10, 1.0, 1], 15, 500, 'gauss_point', 1:500, 3, 0);

%with noise (better)
[cost, rts, w, w_update, u_all]=weightfit([0.10, 1.0, 1], 15, 500, 'gauss_point', 1:500, 3, 0, 15);

parpool('local')
%fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set', 'DiffMinChange', 0.01);
fmincon_options = optimoptions(@fmincon, 'DiffMinChange', 0.01);
fmincon_options.MaxFunEvals = 5000;
fmincon_options.MaxIter = 5000;

lower_bounds=[0, 0.3]; %prop_offset, basis_sep, decay_zero
upper_bounds=[0.5, 2.5];
init_params=[0.125, 1.5];

nbasis=14;
ntrials=1000;
tvec=1:500;

[par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis, ntrials, 'gauss_point', tvec, 3), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

[cost, rts, w, w_update, u_all] = weightfit(par, nbasis, ntrials, 'gauss_point', tvec, 3, 0, 20);
[cost, rts, w, w_update, u_all] = weightfit([.1 2], nbasis, 200, 'gauss_point', tvec, 3, 1, 20);

init_params_1=ones(1,length(centers_1offset));
upper_bounds_1=ones(1,length(centers_1offset)).*100;
lower_bounds_1=ones(1,length(centers_1offset)).*-100;


[par_1offset, cost_1offset, exitflag] = fmincon(@(params) correctrbf(params, 1:500, centers_1offset, ones(1,length(centers_1offset)).*sig, [1 500]), init_params_1, [], [], [], [], lower_bounds_1, upper_bounds_1, [], fmincon_options);
f = rbfeval(1:500, par_1offset, centers_1offset, ones(1,length(centers_1offset)).*sig, [1 500]);


%fmincon_options = optimset(@fmincon);
%fmincon_options.MaxFunEvals = 100000;
%fmincon_options.MaxIter = 100000;

parpool('local')
fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set', 'DiffMinChange', 0.01);
fmincon_options.MaxFunEvals = 5000;
fmincon_options.MaxIter = 5000;

fmincon_options = optimoptions(@fmincon,'DiffMinChange',1);

%test triangle update
lower_bounds=[0, 0.3, 20]; %prop_offset, basis_sep, decay_zero
upper_bounds=[0.5, 2.5, 350];
init_params=[0.125, 1.5, 100];

nbasis=20;
ntrials=1000;
func='triangle';
tvec=1:500;
edge_correct=1; %fix value of temporal generalization at the edge to be the value at the min/max of the observed interval.

%cost and functions at initial values
[cost, rts, w, w_update, u_all] = weightfit(init_params, nbasis, ntrials, func, tvec, edge_correct);


nbasis_test=5:30;
costs_triangle=zeros(length(nbasis_test), 1);
pars_triangle=zeros(length(nbasis_test), length(init_params));

parfor n=1:length(nbasis_test)
    [par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis_test(n), ntrials, func, tvec, edge_correct), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
    costs_triangle(n) = cost;
    pars_triangle(n,:) = par;
end

%print results
horzcat(nbasis_test', costs_triangle,pars_triangle)

%get results at fitted params
[cost, rts, w, w_update, u_all] = weightfit(pars_triangle(5,:), nbasis, ntrials, func, tvec, edge_correct);

%15 basis functions
%15.0000    6.9640    0.1484    1.3957   65.3530
[cost, rts, w, w_update, u_all] = weightfit(pars_triangle(11,:), 15, ntrials, func, tvec, edge_correct);

%24 basis functions with sharp update and highly separarated basis functions... weird
%24.0000    0.9700    0.0014    2.9988   26.0650
[cost, rts, w, w_update, u_all] = weightfit(pars_triangle(20,:), 24, ntrials, func, tvec, edge_correct);

for i = 1:length(nbasis_test)
    [cost, rts, w, w_update, u_all] = weightfit(pars_triangle(i,:), nbasis_test(i), ntrials, func, tvec, edge_correct);
    figure(1); hist(rts);  text(250, 120, sprintf('nbasis: %d', nbasis_test(i)));
    figure(2); plot(u_all(1000,:));
    figure(3); plot(w(1000,:));
    drawnow update;
    pause(0.2);
end

save('tempfunction_costs_pars.mat', 'nbasis_test', 'costs_triangle', 'pars_triangle', 'costs_exp', 'pars_exp');

[cost, rts, w, w_update, u_all] = weightfit(pars_triangle(14,:), nbasis_test(14), ntrials, func, tvec, edge_correct);

%the basis at 24 with no overlap looks bizarre... just capitalizing on chance here and the separation of the basis rel.
%the chi-square bins
figure(1); hist(rts)
figure(2); plot(u_all(1000,:))
figure(3); plot(w(1000,:))



lower_bounds=[0, 0.3, 0.8]; %prop_offset, basis_sep, lambda
upper_bounds=[0.8, 2.5, 0.99999];
init_params=[0.125, 1.5, 0.95];

func='exponential';

%test exponential with 15 basis functions
%[cost, rts, w, w_update, u_all]=weightfit([0.125, 1.5, 0.95], 15, 1000, 'exponential', 1:500, 0);
[cost, rts, w, w_update, u_all]=weightfit([0.125, 1.5, 50], 15, 1000, 'square_auc', 1:500, 0);
[cost, rts, w, w_update, u_all]=weightfit([0.125, 1.2, 0.98], 15, 1000, 'exp_auc', 1:500, 0);

[par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis_test(n), ntrials, func, tvec, 0), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

%exponential update
parfor n=1:length(nbasis_test)
    [par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis_test(n), ntrials, func, tvec, edge_correct), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
    costs_exp(n) = cost;
    pars_exp(n,:) = par;
end

horzcat(nbasis_test', costs_exp',pars_exp)



for i = 1:length(nbasis_test)
    [cost, rts, w, w_update, u_all] = weightfit(pars_exp(i,:), nbasis_test(i), ntrials, func, tvec, edge_correct);
    figure(1); hist(rts); text(250, 120, sprintf('nbasis: %d', nbasis_test(i)));
    figure(2); plot(u_all(1000,:));
    figure(3); plot(w(1000,:));
    drawnow update;
    pause(0.2);
end

%params=[.2 1.0 0.91];
params=[.125 1.5 0.9];
[cost, rts, w, w_update, u_all] = weightfit(params, 20, 150, 'exp_auc', tvec, 1, 1);
hist(rts)

%decent behavior
params=[.125 1.5 50];
[cost, rts, w, w_update, u_all] = weightfit(params, 20, 5000, 'square_auc', tvec, 2, 0);
hist(rts)

%truncated gaussian basis
params=[0 1.5 50];
[cost, rts, w, w_update, u_all] = weightfit(params, 20, 5000, 'square_auc', tvec, 0, 0);
hist(rts)


%square auc update
% lower_bounds=[0, 0.3, 10]; %prop_offset, basis_sep, decay_zero
% upper_bounds=[0.5, 3, 250];
% init_params=[0.125, 1.5, 80];

%fixed separation of 1.5
lower_bounds=[0, 10]; %prop_offset, basis_sep, decay_zero
upper_bounds=[0.5, 250];
init_params=[0.125, 80];

%fixed separation of 1.5, offset of 12.5%
lower_bounds=[10]; %prop_offset, basis_sep, decay_zero
upper_bounds=[250];
init_params=[80];

costs_squareauc=zeros(length(nbasis_test), 1);
pars_squareauc=zeros(length(nbasis_test), length(init_params));
func='square_auc';
parfor n=1:length(nbasis_test)
    [par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis_test(n), 500, func, tvec, edge_correct), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
    costs_squareauc(n) = cost;
    pars_squareauc(n,:) = par;
end

%exp_auc update
params=[.3 1.5 0.99];
[cost, rts, w, w_update, u_all] = weightfit(params, 15, 100, 'exp_auc', tvec, 1, 1);
hist(rts)


%test truncated gauss (no offset)
params=[0 1.5 0.98];
[cost, rts, w, w_update, u_all] = weightfit(params, 15, 100, 'exp_auc', tvec, 0, 1);
hist(rts)



lower_bounds=[0, 0.3, 0.8]; %prop_offset, basis_sep, lambda
upper_bounds=[0.5, 2.5, 0.9999];
init_params=[0.125, 1.5, 0.95];

fmincon_options = optimoptions(@fmincon, 'DiffMinChange', 0.01);
func='expauc';
nbasis_test=10:25;
costs_expauc=zeros(length(nbasis_test), 1);
pars_expauc=zeros(length(nbasis_test), length(init_params));
func='square_auc';
parfor n=1:length(nbasis_test)
    [par, cost, exitflag] = fmincon(@(params) weightfit(params, nbasis_test(n), 500, func, tvec, edge_correct), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
    costs_squareauc(n) = cost;
    pars_squareauc(n,:) = par;
end

mov=repmat(struct('cdata', [], 'colormap', []), size(u_all, 1),1);
for i=1:size(u_all,1)
    plot(u_all(i,:), 'r'); hold on;
    if i > 1, plot(u_all(i-1,:), 'b'); end
    scatter(rts(i), 0.2*max(u_all(i,:)), 200);
    hold off;
    drawnow update;
    
    mov(i) = getframe(gcf);
end

movie2avi(mov, '~/test.avi', 'compression', 'None', 'FPS', 5);
    

movie = immovie(X,map); #% map is the colormap you want to use

implay(movie);

%look at effect of offset
% offsets=0:0.01:1.0;
% costs=[];
% params=[]
% for o=1:length(offsets)
%     [par, cost, exitflag] = fmincon(@(params) rbffit(params, x, y, 1, sds, offsets(o)), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
%     costs(o) = cost;
%     params = vertcat(params, par);
% end
% 
% %re-run at fitted values to get centers
% %use best offset from above
% offset_best=offsets(find(costs == min(costs)));
% fprintf('best offset: %.3f\n', offset_best);
% params_best=params(find(costs == min(costs)), :);
% [cost, yhat, centers] = rbffit(params_best, x, y, 1, sds, offset_best);



%look at effect of sd
sds=0.05:.01:1.0;
offsets=0.5; %fix at reasonable value
costs=[];
params=[]
for s=1:length(sds)
    [par, cost, exitflag] = fmincon(@(params) rbffit(params, x, y, 1, sds(s), offset), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
    costs(s) = cost;
    params = vertcat(params, par);
end

%re-run at fitted values using best params
sd_best=sds(find(costs==min(costs)));
fprintf('best sd: %.3f\n', sd_best);
params_best=params(find(costs == min(costs)), :);
[cost, yhat, centers] = rbffit(params_best, x, y, 1, sd_best, offset);

%extract fitted parameters
sd_fitted=sd_best;
centers_fitted=centers;
weights_fitted=params_best;

%model 2
%sd_fitted = params(1);
%weights_fitted=params(2:end);

%model 3
%sd_fitted = params(1);
%centers_fitted=params(2:((length(params)-1)/2)+1);
%weights_fitted=params(((length(params)-1)/2)+2:end);


figure(1)
scatter(x,y)
hold on
scatter(x,yhat, 'r')

%plot rbf functions

%construct radial basis matrix using Gaussians
t=500;
gaussmat = zeros(nbasis,500);
xrng=xmin:(xmax-xmin)/(t-1):xmax; %x axis for plotting 
for j = 1:nbasis
    gaussmat(j,:) = gaussmf(xrng,[sd_fitted centers_fitted(j)]);
end

figure(2);
rbf_mat = weights_fitted'*ones(1,t) .* gaussmat;
plot(xrng,rbf_mat);
ylabel('temporal basis function')


rbf_sum = rbfeval(xrng, weights_fitted, centers_fitted, sd_fitted);
figure(3);
plot(xrng, rbf_sum)

%optim_options = optimset(@fminsearch);
%[fitted, cost, yhat] = fminsearch(@(params) rbffit(params, x, y), params, optim_options);

%plot fitted function
%pred=rbfeval(x, fitted(2:end));
