function diagnose_optimality(model)
%cd '/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing/output';
cd '/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/optimality_testing';

%model='fixedLR_decay';
%model='kalman_uv_sum';
%model='fixedLR_softmax';
optfiles = glob(['output/optimize_output_sinusoid_', model, '*.mat']);
disp(optfiles);

nreps = 10;
nbest = 10; %number of parameter sets to test
ntimesteps=500;
ntrials = 60;
cliffpullback=20;

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

rng(625); %fix seed for pulling reward probabilities
%randomly sample with replacement the desired number of contingencies
keep = randsample(1:ntimesteps, nreps, true);

clear optmat;
for k = 1:length(keep)
    thisCont=[];
    thisCont.name = ['sinusoid' num2str(keep(k))];
    thisCont.sample = zeros(1, ntimesteps); %keeps track of how many times a timestep has been sampled by agent
    thisCont.lookup = zeros(ntimesteps, ntrials); %lookup table of timesteps and outcomes
    thisCont.ev = allshift(keep(k),:,1);
    thisCont.prb = allshift(keep(k),:,2);
    thisCont.mag = allshift(keep(k),:,3);
    
    rvec = rand(ntimesteps, ntrials);
    for t = 1:ntrials
        thisCont.lookup(:,t) = (allshift(keep(k),:,2) > rvec(:,t)') .* allshift(keep(k),:,3);
    end

    thisCont.seeds = randi([1 500], 1, 5);
    thisCont.firstrt = randi([1 500], 1, 1);

    optmat(k) = thisCont;
end

nagents=length(optfiles);

costs=NaN(nagents, nreps, nbest);
%v_it=NaN(nagents, nreps, ntrials, ntimesteps);
rts = NaN(nagents, nreps, nbest, ntrials);
%ret = cell(nagents, nreps); 

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

for o = 1:(length(optfiles) + 1)
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
    
    rng(1040); %consistent starting point for randomization of contingency

    %pars = cell2mat(opt.pars);
    %medpars = median(pars);
    %bestpars = pars(find(opt.costs==min(opt.costs)),:);

    for q = 1:nbest
        parsq = candidates(q, 1:(size(candidates,2)-1)); %trim last column, which contains cost
        costs=NaN(nreps,1);
        ret=cell(nreps,1);
        for p = 1:nreps
            if strcmpi(a.name, 'franktc')
                [costs(p), ~, ret{p}] = TC_Alg_forward(parsq, priors, optmat(p), optmat(p).seeds, a.ntrials, [0 5000]);
            else
                %parsq=[.04, 1.8, .08, .2];
                [costs(p), ~, ~, ret{p}] = clock_sceptic_agent(parsq, a, optmat(p).seeds, optmat(p), a.ntrials, a.nbasis, a.ntimesteps, 20);
                %underlying contingency versus choices (stars)
            end
        end
        sceptic_movie(costs, a, ret, candidates(1:nbest, 1:(size(candidates,2)-1)), optmat, sprintf('movies/%s_pars%d', model, q));
    end
end

end
