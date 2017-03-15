%output relevant stats from vba fits for fmri

function [res] = extract_sceptic_fmri(results_dir)

if nargin < 1, results_dir = '/Users/michael/ics/temporal_instrumental_agent/clock_task/vba_fmri/vba_out'; end

[~, model_description] = fileparts(results_dir);

addpath('/storage/group/mnh5174_collab/temporal_instrumental_agent'); %has glob.m
fitfiles = glob([results_dir, '/*.mat']);
nruns=8;
ntrials = 50; %number of trials in a run (used for dividing according to multisession)
nbasis = 24; %16;
nstates = 3; %V, PE, decay
multisession=0; %multisession requires one to look at the right basis position within a given run (run-related offsets)
npars = 4; %alpha, gamma, prop_spread, beta
  
nsubjs = length(fitfiles);

%runs x trials x nbasis
V = NaN(nsubjs, nruns, nbasis, ntrials);
PE = NaN(nsubjs, nruns, nbasis, ntrials);
D = NaN(nsubjs, nruns, nbasis, ntrials);
pars = NaN(nsubjs, npars+2);

for f = 1:length(fitfiles)
    load(fitfiles{f});
    nruns=size(posterior.muX, 2)/ntrials;
    %for multisession, the muX matrix is 384 x 400
    %400 is the number of trials
    %384 is 16 (basis) x 3 (V, PE, decay) x 8 runs
    for run = 1:nruns
        if multisession    
            runoffset = (run-1)*nbasis*nstates + 1;
        else
            runoffset = 1; %no need to offset the hidden state vector by run
        end
        
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

    [~, fname] = fileparts(fitfiles{f});
    id = str2double(fname(isstrprop(fname,'digit')));
    pars(f,1:(length(posterior.muTheta) + length(posterior.muPhi) + 2)) = [id, out.F, posterior.muTheta', posterior.muPhi];
end

res.fitfiles = fitfiles;
res.V = V;
res.PE = PE;
res.D = D;
res.pars=pars;

save(sprintf('posterior_states_decay_nomultisession_%s.mat', model_description), 'fitfiles', 'V', 'PE', 'D', 'pars');
