cd '/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing/output';

uvsum = load('optimize_output_kalman_uv_sum.mat');

pars = median(cell2mat(uvsum.pars(:,1,6))); %median of each parameter over optimization values

IEVmaster = uvsum.optmat{2}(1); %first contingency is IEV

a = uvsum.agents; %agent used for these parameters

rng(1040);
nreps = 100;
costs=NaN(nreps, 1);
v_it=NaN(nreps, a.ntrials, a.ntimesteps);
rts = NaN(nreps, a.ntrials);

for p = 1:nreps
    seeds = randi([1 500], 1, 4);
    permIEV = IEVmaster;
    permIEV.lookup = permIEV.lookup(:, randperm(size(permIEV.lookup,2))); %randomly permute columns
    
    [costs(p), v_it(p,:,:), rts(p,:), ret(p)] = clock_sceptic_agent(pars, a.name, seeds, permIEV, a.ntrials, a.nbasis, a.ntimesteps);

end

plot(rts')
histogram(mean(rts'))

%%

frank = load('optimize_output_franktc.mat');

pars = median(cell2mat(frank.pars(:,1,6))); %median of each parameter over optimization values

IEVmaster = frank.optmat{2}(1); %first contingency is IEV

a = frank.agents; %agent used for these parameters

rng(1040);
nreps = 100;
costs=NaN(nreps, 1);
v_it=NaN(nreps, a.ntrials, a.ntimesteps);
rts = NaN(nreps, a.ntrials);
priors.V = 0; %don't give agent any insight into previous values (no SCEPTIC receives this boost)
priors.Go = 0;
priors.NoGo = 0;

for p = 1:nreps
    seeds = randi([1 500], 1, 4);
    permIEV = IEVmaster;
    permIEV.lookup = permIEV.lookup(:, randperm(size(permIEV.lookup,2))); %randomly permute columns
    
    [costs(p), rts(p,:), ret(p)] = TC_Alg_forward(pars, priors, permIEV, seeds, a.ntrials, [0 5000]);

end

plot(rts')
histogram(mean(rts'))


%%GRW approach to random contingencies
%%Magnitude varies between 100 and 500
mag_start = 300;
step_size_mag = 100;
mag_max=500;
mag_min=100;

%%probability varies between 0.3 and 0.7
step_size_prb = 0.15;
prb_start = 0.5;
prb_max=0.7;
prb_min=0.3;

for rep=1:100

ntimesteps=500;

steps_mag=randn(ntimesteps,1)*step_size_mag;
steps_prb=randn(ntimesteps,1)*step_size_prb;

mag = zeros(ntimesteps,1);
prb = zeros(ntimesteps,1);
mag(1) = mag_start;
prb(1) = prb_start;
for i = 2:ntimesteps
    mag(i) = mag(i-1) + steps_mag(i);
    prb(i) = prb(i-1) + steps_prb(i);
end

%normalize to desired range
mag_scaled = (mag - min(mag))*(mag_max-mag_min)/(max(mag)-min(mag)) + mag_min;
prb_scaled = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min;

figure(1); clf; subplot(3,1,1); plot(smooth(mag_scaled, 40));
subplot(3,1,2); plot(smooth(prb_scaled, 30));
subplot(3,1,3); plot(smooth(prb_scaled.*mag_scaled, 40));

drawnow update;
pause(1);

end






%plot(smooth(rand(500, 1), 100))


