%test_rbf script

xy = dlmread('xy_cor0.9.txt', '\t');
x=xy(:,1);
y=xy(:,2);
nbasis=12;

%perfect straight line
%best offset: 0.49, sd = 0.3, nbasis = 12
%best sd: 0.65, offset = 0.5, nbasis = 12 (but very flat costs for sd > 0.29)
%x=[-2:.1:2]';
%y=[-2:.1:2]';

%perfect quadratic
%best offset: 0.53
%best sd: 0.31, offset = 0.5, nbasis = 12
%x=[-2:.1:2]';
%y=x.^2;

%perfect cubic
%best offset: 0.58, sd = 0.3, nbasis = 12
%best sd: 0.31, offset = 0.5, nbasis = 12
%x=[-2:.1:2]';
%y=x.^3;

%perfect quartic
%best offset: 0.63, sd = 0.3, nbasis = 12
%best sd: 0.31, offset = 0.5, nbasis = 12
%x=[-2:.1:2]';
%y=x.^4;

%sine function
%best offset: 0.08, sd = 0.3
%best sd: 0.59, offset = 0.5, nbasis=12 
%the function has about 3 cycles here, so to match the complexity, a
%smaller offset helps with giving enough basis functions inside the
%observed interval. If we increase the number of basis functions to 16 or
%20, the costs become pretty flat (i.e., most offsets do well enough.
x=[-2:.05:2]';
y=sin(x*5)*3;

%conclusions about offset: offsets around 0.5 (12.5% of the interval) on each side are
%optimal for placing centers. This gives enough room at the edges to fit
%almost any function. If the function becomes highly complex, more basis
%functions will be needed and/or the offset must be decreased to give
%more accurate support for the function inside the observed interval.

%conclusion about sd: for polynomials up to x^4, sd around 0.30 gives excellent results. This is in the
%context of 12 basis functions spread across the interval -2.5 -- 2.5. Thus, basis functions are spaced 0.45
%units apart. With an SD of 0.3, the Gaussians are spread Cohen's d = 0.45 / 0.3 = 1.52. A Cohen's d of 1.52
%corresponds to 45% overlap between distributions. Thus, having about 50% overlap between basis functions
%seems optimal for approximating most functions.

fmincon_options = optimset(@fmincon);
fmincon_options.MaxFunEvals = 100000;
fmincon_options.MaxIter = 100000;

%calculate centers of basis functions on the basis of x (time) scaling
xmin=min(x);
xmax=max(x);


% lower_bounds=[0 (xmin-1)*ones(1,nbasis) -100*ones(1,nbasis)];
% upper_bounds=[100 (xmax+1)*ones(1,nbasis) ones(1,nbasis)*100];
% init_params=[0.35 centers zeros(1,nbasis)]; %initial parameters

% lower_bounds=[0 -100*ones(1,nbasis)];
% upper_bounds=[100 ones(1,nbasis)*100];
% init_params=[0.35 zeros(1,nbasis)]; %initial parameters

%weights only (model 1)
lower_bounds=[-100*ones(1,nbasis)];
upper_bounds=[ones(1,nbasis)*100];
init_params=[zeros(1,nbasis)]; %initial parameters

sds=0.3; %fixed for model 1
offset=0.5; %space rbfs a bit beyond observed data
model=1; %sds fixed, centers fixed, weights free
centers=(xmin-offset):(xmax-xmin)/(nbasis-1):(xmax+offset);


[cost, yhat]=rbffit(init_params, x, y); %fit at initial settings

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
