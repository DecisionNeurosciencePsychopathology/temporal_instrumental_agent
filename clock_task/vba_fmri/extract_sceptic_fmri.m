%output relevant stats from vba fits for fmri
results_dir = '/Users/michael/vba_out'; %output directory
fitfiles = glob([results_dir, '/*.mat']);
nruns=8;
ntrials = 50; %number of trials in a run (used for dividing according to multisession)
nbasis = 16;
nstates = 3; %V, PE, decay

nsubjs = length(fitfiles);

%runs x trials x nbasis
V = NaN(nsubjs, nruns, nbasis, ntrials);
PE = NaN(nsubjs, nruns, nbasis, ntrials);
D = NaN(nsubjs, nruns, nbasis, ntrials);

for f = 1:length(fitfiles)
    load(fitfiles{f});
    %the muX matrix is 384 x 400
    %400 is the number of trials
    %384 is 16 (basis) x 3 (V, PE, decay) x 8 runs
    for run = 1:nruns
        runoffset = (run-1)*nbasis*nstates + 1;
        trialoffset = (run-1)*ntrials;
        for state = 1:nstates
            stateoffset = (state-1)*nbasis;
            muXslice = posterior.muX((runoffset+stateoffset):(runoffset+stateoffset+nbasis - 1), (trialoffset+1):(trialoffset+ntrials));
            if state == 1 %value
                V(f, run, :, :) = muXslice;
            elseif state == 2 %prediction error
                PE(f, run, :, :) = muXslice;
            elseif state == 3 %decay
                D(f, run, :, :) = muXslice;
            end
        end
    end
end

save('posterior_states_decay.mat', 'fitfiles', 'V', 'PE', 'D');