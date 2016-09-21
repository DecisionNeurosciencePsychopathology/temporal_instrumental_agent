function compute_basis_adjustment(ntimesteps, nbasis, prop_spread, decay, peversion, vchoice)
addpath('../')
if nargin < 4, decay = 0; end
if nargin < 5, peversion = 0; end %regular PE update alpha*e(reward - value)
if nargin < 6, vchoice = 0; end %value-guided choice (1) or random explore (0)

%use a fixed learning rate model with the basis to figure out how to reach asymptote for each basis function
%what adjustments have to be made to the basis functions to get there?
[~, ~, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(ntimesteps, nbasis, prop_spread);
alpha=.05;
beta = 10; %typical softmax inverse temperature
%draw a random sample from the phase-shifted sinusoids for each replication
ev = 10*sin(2*pi*(1:ntimesteps).*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps);
ev = ev + abs(min(ev)) + 10;
prb = 25*cos(2*pi*(1:ntimesteps).*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps);
prb_max=0.7;
prb_min=0.3;
prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min;
mag = ev./prb;

%shift sinusoid randomly each call
offset=randi([0, ntimesteps], 1, 1);

shift=[offset:ntimesteps 1:(offset-1)];
ev = ev(shift);
prb = prb(shift);
mag = ev./prb;

%for now, always deliver reward (no probability)
tolerance = 2; %2 point discrepancy between basis approximation and underlying contingency
mu_ij = zeros(nbasis,1);
v_func = zeros(1, ntimesteps);

while(any(abs(v_func - mag) > tolerance))
    %sample each timestep in random order (so that each is selected)
    
    if vchoice == 0
        pvec=randperm(ntimesteps);
        for t = 1:ntimesteps
            choice = pvec(t);
            mu_ij = updatev(choice, mu_ij);
            plotv(choice); %also recomputes v_func
        end
        
        interimplot;
    else
        %value max choice: Fixed LR V, potentially with decay
        %arbitrary 100 trials
        choice = randi(ntimesteps, 1, 1); %random initial choice
        for i = 1:100
            mu_ij = updatev(choice, mu_ij);
            plotv(choice); %also recomputes v_func
            
            %softmax choice rule over value function
            pchoice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature
            %if (all(v_func==0)), v_func=rand(1, length(v_func)).*1e-6; end; %need small non-zero values to unstick softmax on first trial
            choice = randsample(tvec, 1, true, pchoice); %softmax_stream,
        end
        
        interimplot;
    end
    
end

%%Nested functions below for plotting and value update
    function plotv(choice)
        figure(1);
        plot(1:ntimesteps, mag, 'g');
        hold on;
        plot(1:ntimesteps, v_func, 'b'); %t-1
        v_jt=mu_ij*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
        v_func = sum(v_jt); %subjective value by timestep as a sum of all basis functions
        
        plot(1:ntimesteps, v_func, 'r');
        plot(choice, mag(choice) + 1,'*b');
        hold off
        drawnow update;
    end

    function interimplot
        %after a set of trials, display correlation on plot and pause
        r = corr(v_func', mag')
        text(ntimesteps - 50, max([v_func, mag]) - 10, sprintf('%.3f', r));
        
        %also show a plot of the difference in a separate panel
        figure(2); 
        dvec = v_func - mag;
        [~, msort] = sort(mag);
        rdiff = corr(mag', v_func' - mag')
        scatter(mag(msort), dvec(msort));
        
        %compute reference line (apparently the intercept of 0 is actually the min in MATLAB)
        %refline([range(dvec)/range(mag), min(dvec)]);
        refline([range(dvec)/range(mag), 0]);
        text(max(mag) - 15, max(v_func - mag) - 5, sprintf('r(mag, diff) = %.3f', rdiff));
        %text(ntimesteps - 150, max(v_func - mag) - 10, sprintf('r(mag, diff) = %.3f', rdiff));
        pause(1.2);
    end

    function newmu = updatev(choice, curmu)
        %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
        elig = gaussmf(tvec, [sig_spread, choice]);
        
        %compute sum of area under the curve of the gaussian function
        auc=sum(elig);
        elig=elig/auc*refspread;
        e_ij = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 2);
        
        if peversion == 0
            %PE version 1: familiar e*(rew - Vb)
            delta_ij = e_ij.*(mag(choice) - curmu);
        elseif peversion == 1
            %PE version 2: based on (our understanding of) discussion with Yael Niv: rew - e*Vfunc
            
            %estimated value with respect to time
            %sumv_ij = sum(curmu*ones(1,ntimesteps) .* gaussmat, 2);
            
            %compute value as a function of time (evaluating basis)
            v_i=sum(curmu*ones(1,ntimesteps) .* gaussmat);
            v_it = v_i(choice); %scalar value estimate at chosen response
            
            %now we want the prediction error for each basis (WRONG)
            %delta_ij = mag(pvec(t)) - ones(nbasis,1)*elig .* v_ij;
            
            %distribute estimated value at RT according to eligibility (works correctly)
            delta_ij = e_ij.*(mag(choice) - v_it);
        end
        
        d = -decay.*(1-e_ij).*curmu;
        newmu = curmu + alpha.*delta_ij + d;
    end

end

