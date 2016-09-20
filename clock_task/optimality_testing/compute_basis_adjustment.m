function compute_basis_adjustment(ntimesteps, nbasis, prop_spread, decay, peversion)
addpath('../')
if nargin < 4, decay = 0; end
if nargin < 5, peversion = 0; end %regular PE update alpha*e(reward - value)

%use a fixed learning rate model with the basis to figure out how to reach asymptote for each basis function
%what adjustments have to be made to the basis functions to get there?
[~, ~, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(ntimesteps, nbasis, prop_spread);
alpha=.05;
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
tolvec = zeros(1,ntimesteps);
mu_ij = zeros(nbasis,1);
v_func = zeros(ntimesteps,1);
while(any(abs(tolvec - mag) > tolerance))
    %sample each timestep in random order
    pvec=randperm(ntimesteps);
    for t = 1:ntimesteps
        %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
        elig = gaussmf(tvec, [sig_spread, pvec(t)]);
        
        %compute sum of area under the curve of the gaussian function
        auc=sum(elig);
        elig=elig/auc*refspread;
        e_ij = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 2);
        
        
        if peversion == 0
            %PE version 1: familiar e*(rew - Vb)
            delta_ij = e_ij.*(mag(pvec(t)) - mu_ij);
        elseif peversion == 1
            %PE version 2: based on (our unnderstanding of) discussion with Yael Niv: rew - e*Vfunc
            
            %estimated value with respect to time
            %sumv_ij = sum(mu_ij*ones(1,ntimesteps) .* gaussmat, 2);
            
            v_i=sum(mu_ij*ones(1,ntimesteps) .* gaussmat); %use vector outer product to replicate weight vector
            v_it = v_i(pvec(t)); %scalar value estimate at chosen response
            
            %now we want the prediction error for each basis (WRONG)
            %delta_ij = mag(pvec(t)) - ones(nbasis,1)*elig .* v_ij;
            
            %distribute estimated value at RT according to eligibility (works correctly)
            delta_ij = e_ij.*(mag(pvec(t)) - v_it);
        end
        
        
        d = -decay.*(1-e_ij).*mu_ij;
        mu_ij = mu_ij + alpha.*delta_ij + d;
        
        plot(1:ntimesteps, mag, 'g');
        hold on;
        plot(1:ntimesteps, v_func, 'b'); %t-1
        v_jt=mu_ij*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
        v_func = sum(v_jt); %subjective value by timestep as a sum of all basis functions
        
        plot(1:ntimesteps, v_func, 'r');
        plot(pvec(t), mag(pvec(t)) + 1,'*b');
        hold off
        drawnow update;
    end
    
    r = corr(v_func', mag')
    text(ntimesteps - 50, max([v_func, mag]) - 10, sprintf('%.3f', r));
    pause(1);
    

end
    

end
