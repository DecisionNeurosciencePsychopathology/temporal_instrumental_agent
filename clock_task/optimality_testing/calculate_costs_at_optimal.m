%cd '/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing/output';
cd '/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/optimality_testing';
addpath('../');
%optfiles = glob('output/optimize_output*.mat');
optfiles = glob('output/optimize_output_allequate*.mat');

nreps = 100;
nbest = 5; %number of parameter sets to test
ntimesteps=500;
ntrials = [30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
cliffpullback=20;
plots=false;
target='allequate';

if strcmpi(target, 'sinusoid')
    %draw a random sample from the phase-shifted sinusoids for each replication
    ev = 10*sin(2*pi*(1:ntimesteps).*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps);
    ev = ev + abs(min(ev)) + 10;
    prb = 25*cos(2*pi*(1:ntimesteps).*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps);
    prb_max=0.7;
    prb_min=0.3;
    prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min;
    
    allshift = NaN(ntimesteps, ntimesteps, 3);
    conds = 1:ntimesteps;
    for i = 1:ntimesteps
        shift=[i:ntimesteps 1:(i-1)];
        evi = ev(shift);
        prbi = prb(shift);
        
        allshift(i,:,1) = evi;
        allshift(i,:,2) = prbi;
        allshift(i,:,3) = evi./prb;
    end
    
    %randomly sample with replacement the desired number of contingencies
    keep = randsample(1:ntimesteps, nreps, true);

    rng(625); %fix seed for pulling reward probabilities
    
    clear optmat;
    for k = 1:length(keep)
        thisCont=[];
        thisCont.name = ['sinusoid' num2str(keep(k))];
        thisCont.sample = zeros(1, ntimesteps); %keeps track of how many times a timestep has been sampled by agent
        thisCont.lookup = zeros(ntimesteps, max(ntrials)); %lookup table of timesteps and outcomes
        thisCont.ev = allshift(keep(k),:,1);
        thisCont.prb = allshift(keep(k),:,2);
        thisCont.mag = allshift(keep(k),:,3);
        
        rvec = rand(ntimesteps, max(ntrials));
        for t = 1:max(ntrials)
            thisCont.lookup(:,t) = (allshift(keep(k),:,2) > rvec(:,t)') .* allshift(keep(k),:,3);
        end
        
        optmat(k) = thisCont;
    end

elseif strcmpi(target, 'allequate')
    %Aug2016: four core contingencies (IEV, DEV, CEV, CEVR) with equal EV AUCs (to weigh equally in optimization)
    load('mastersamp_equateauc.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
    condnames={'IEV', 'DEV', 'CEV', 'CEVR'}; %by default, optimize over all contingencies in lookup
    ncond=length(condnames); %how many conditions/contingencies
    ntimesteps = 500;
    
    repspercond = floor(nreps/length(condnames));
    condmat = repmat(condnames, repspercond, 1);
    nrep = cell(1, nreps);
    for i = 1:ncond
        for j = 1:repspercond
            nrep{(i-1)*repspercond + j} = num2str(j);
        end
    end
    %condnames = strcat(condmat(:)', nrep);
    condrep = condmat(:)';
    
    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    clear optmat;
    for i = 1:ncond
        for j = 1:repspercond
            tmp = mastersamp_equateauc.(condnames{i});
            tmp.name = condnames{i}; %for keeping track of what contingency this is (esp. when they mix in ALL)
            tmp.lookup = tmp.lookup(1:ntimesteps, randperm(size(tmp.lookup,2))); %randomly permute columns
            tmp.ev = tmp.ev(1:ntimesteps);
            tmp.sample = tmp.sample(1:ntimesteps);
            optmat((i-1)*repspercond + j) = tmp; %add to a vector of structs
        end
    end
end

nagents=length(optfiles);

priors.V = 0; %don't give agent any insight into previous values (no SCEPTIC receives this boost)
priors.Go = 0;
priors.NoGo = 0;

%fly through the contingencies (looks good)
%for i = 1:length(optmat)
%    pause(0.08);
%    subplot(3,1,1);
%    plot(optmat(i).prb);
%    subplot(3,1,2);
%    plot(optmat(i).mag);
%    subplot(3,1,3);
%    plot(optmat(i).ev);
%end

%storage structure for output
oresults=struct();

for o = 1:(length(optfiles) + 1)
    if o <= length(optfiles);
        opt = load(optfiles{o});
        a=opt.agents;
        
        %look at top 10% and compare parameters
        bestqcosts = sort(opt.costs);
        bestqcosts = bestqcosts(1:nbest);
        %cutoff = quantile(opt.costs, .1); %bottom 10% since lower costs are better
        cutoff = max(bestqcosts);
        bestcost_indices = opt.costs <= cutoff;
        bestcosts = opt.costs(bestcost_indices);
        %[~, ind] = sort(bestcosts);
        candidates = cell2mat(opt.pars(bestcost_indices));
        candidates = horzcat(candidates, bestcosts);
        candidates = sortrows(candidates, size(candidates, 2));
        %candidates = candidates(ind, :); %sort in terms of low to high (good to bad) costs)
        
        %plot parameters for best 10 optimizations
        if plots
            fig = figure;
            npars = size(candidates, 2);
            for jj = 1:npars
                subplot(npars, 1, jj);
                if jj==npars
                    plot(1:size(candidates, 1), candidates(:,jj));
                    title('total cost');
                else
                    bar(candidates(:,jj));
                    title(sprintf('%s',a.parnames{jj}), 'Interpreter', 'none');
                end
            end
            [~, fbase, ~] = fileparts(optfiles{o});
            
            print(fig,['output/' fbase '_best10.png'],'-dpng')
            
            %plot histograms of parameters for all optimizations
            fig = figure;
            allpars = cell2mat(opt.pars);
            allpars = horzcat(allpars, opt.costs);
            allpars = sortrows(allpars, size(allpars, 2));
            
            npars = size(allpars, 2);
            for jj = 1:npars
                subplot(npars, 1, jj);
                if jj==npars
                    plot(1:size(allpars, 1), allpars(:,jj));
                    title('total cost');
                else
                    hist(allpars(:,jj));
                    title(sprintf('%s',a.parnames{jj}), 'Interpreter', 'none');
                end
            end
            [~, fbase, ~] = fileparts(optfiles{o});
            
            print(fig,['output/' fbase '_allopt.png'],'-dpng')
        end
    else
        a.name = 'null';
        a.nbasis = 24;
        a.ntimesteps = 500;
        a.prop_spread=0; %need this to avoid problem with parameter loading
        a.fixps=1; a.fixbeta=0; %hack: fix prop spread to have the parameter lookup work
        parsq=[];
    end
        
    rng(1040); %consistent starting point for randomization of contingency

    %pars = cell2mat(opt.pars);
    %medpars = median(pars);
    %bestpars = pars(find(opt.costs==min(opt.costs)),:);

    %v_it=NaN(nagents, nreps, ntrials, ntimesteps);
    %rts = NaN(length(ntrials), nreps, nbest, ntrials);
    %ret = cell(nagents, nreps);

    oresults.(a.name).costs = NaN(length(ntrials), nreps, nbest);
    oresults.(a.name).pars  = cell(length(ntrials), nreps, nbest);
    oresults.(a.name).rts   = cell(length(ntrials), nreps, nbest);
    oresults.(a.name).vcorr = cell(length(ntrials), nreps, nbest);
    for t = 1:length(ntrials)
        %rts   = NaN(ntrials(t),1);
        vcorr = NaN(ntrials(t),1);

        for p = 1:nreps
            seeds = randi([1 500], 1, 5);
            optmat(p).firstrt = randi([1 500], 1, 1);
            for q = 1:nbest
                parsq = candidates(q, 1:(size(candidates,2)-1)); %trim last column, which contains cost
                oresults.(a.name).pars{t,p,q} = parsq;
                if strcmpi(a.name, 'franktc')
                    %[costs(o,p), rts(o,p,:), ret{o,p}] = TC_Alg_forward(bestpars, priors, optmat(p), seeds, a.ntrials, [0 5000]);
                    [oresults.(a.name).costs(t,p,q), oresults.(a.name).rts{t,p,q}, ~] = TC_Alg_forward(parsq, priors, optmat(p), seeds, ntrials(t), [0 5000]);
                elseif strcmpi(a.name, 'qlearning')
                    a.clock_options.episodeCount = ntrials(t);
                    a.clock_options.ntimesteps = a.ntimesteps/10; %TD models operate on 50 timesteps, not 500, typically
                    [oresults.(a.name).costs(t,p,q), ~, oresults.(a.name).rts{t,p,q}] = ClockWalking_3D_discountedEv_optimize(a.clock_options, optmat(p), seeds, parsq);
                else
                    %[costs(o,p), v_it(o,p,:,:), rts(o,p,:), ret{o,p}] = clock_sceptic_agent(bestpars, a, seeds, optmat(p), ntrials, a.nbasis, a.ntimesteps);
                    %[costs(o,p), ~, rts(o,p,:), ret{o,p}] = clock_sceptic_agent(bestpars, a, seeds, optmat(p), ntrials, a.nbasis, a.ntimesteps);
                    [oresults.(a.name).costs(t,p,q), ~, oresults.(a.name).rts{t,p,q}, ret] = clock_sceptic_agent(parsq, a, seeds, optmat(p), ntrials(t), a.nbasis, a.ntimesteps, cliffpullback);
                    
                    for kk=1:size(ret.v_it,1)
                        vcorr(kk) = corr(ret.v_it(kk,:)', optmat(p).ev');
                        if plots
                            plot(ret.v_it(kk,:)');
                            hold on
                            plot(optmat(p).ev', 'r');
                            title(sprintf('model: %s, i:%d, rt_i: %d', a.name, kk, ret.rts(kk)));
                            pause(0.1);
                            hold off
                        end
                    end
                    oresults.(a.name).vcorr{t,p,q} = vcorr;
                end
            end            
        end
    end
end

save('/Users/michael/performance_at_optimal.mat', 'oresults', 'optmat', 'optfiles'); %'-v7.3'
%this makes the file about 8GB!  'ret' 'v_it', 
%save('/Users/michael/bootstrap_costs_at_best10.mat', 'costs', 'rts', 'optmat', 'optfiles'); %'-v7.3'

return




%old stuff
%uvsum = load('optimize_output_kalman_uv_sum.mat');

%uvsum = load('optimize_output_sinusoidsingle_kalman_uv_sum');
%uvsum = load('optimize_output_sinusoidsingle_fixedLR_softmax_bfix.mat');
uvsum = load('optimize_output_sinusoidsingle_fixedLR_softmax_1run_freebeta.mat');
%pars = median(cell2mat(uvsum.pars(:,1,6))); %median of each parameter over optimization values
pars = median(cell2mat(uvsum.pars)); %median of each parameter over optimization values

%IEVmaster = uvsum.optmat{2}(1); %first contingency is IEV
IEVmaster = uvsum.optmat{1}(1); %first contingency is IEV

a = uvsum.agents; %agent used for these parameters

rng(1040);
nreps = 100;
costs=NaN(nreps, 1);
v_it=NaN(nreps, a.ntrials, a.ntimesteps);
rts = NaN(nreps, a.ntrials);

parmat = cell2mat(uvsum.pars);
hist(parmat(:,1)) %prop spread
hist(parmat(:,2)) %beta
hist(parmat(:,3)) %alpha


hist(uvsum.costs)

%pars=[.01, .01, .95];


plot(IEVmaster.ev)
plot(rts')
for i=1:100
    plot(rts(i,:))
    axis([1 100 1 500]);
    pause(0.2)
end

plot(rts(10,:))
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


