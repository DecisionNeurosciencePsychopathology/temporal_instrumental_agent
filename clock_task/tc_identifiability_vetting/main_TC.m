%clear all

%initialize global definitions
globdefs

%initialize behavior of fitting function
%model = 'noemo'; %no parameters vary by emotion
%model = 'emoexplore'; %explore epsilon varies by emotion
%model = 'emogonogo'; %go and no go learning rates vary by emotion

generative = 0; % generative model that makes its own choices and gets reward rather than fitting subject data.
multstart = 1; % use multiple starting points for gradient descent.

%initialize optimizer settings
options = optimset(@fmincon);
%options = optimset(options, 'LargeScale', 'off');

%get a list of subject data to fit

%subject directory is relative to this fit_behavior directory.
subjdir='../subjects/';
subjfiles=dir(strcat(subjdir, '*.mat'));

if isempty(subjfiles)
    error(strcat('Could not find any subject data files in: ', subjdir));
end

subjdata=cell(length(subjfiles),1);
subjconcat=[]; %stores all subject data in one array
subjids=[];
for f = 1:length(subjfiles)
    [subjdata{f}, subjids(f)] = loadTCRun(strcat(subjdir, subjfiles(f).name));
    subjconcat = vertcat(subjconcat, subjdata{f});  
end

refit=0; %whether to refit a subject who has already been run.

%loop over and fit all models
models = { 'noemo', 'emoexplore', 'noemosticky' }; %, 'emogonogo', 'emonogo', ...
    %'noemosticky', 'emoexploresticky', 'emogonogosticky', 'emonogosticky', ...
    %'noemo_scram', 'noemosticky_scram' };

%matlabpool(4)

for m = 1:length(models)
    model = models{m};
    bestFit_all = [];
    
    %for now, just fit per subject
    for f = 1:size(subjdata,1)
    %parfor f = 1:size(subjdata,1)
        %skip
        if (exist(strcat('../outputs/parameter_mat/modelVars_', num2str(subjids(f)), model, '.mat'), 'file') > 0 && ~refit)
            fprintf('Skipping: modelVars_%d%s.mat\n', subjids(f), model);
            load(strcat('../outputs/parameter_mat/modelVars_', num2str(subjids(f)), model));
            bestFit_all(f, :) = bestFit;
            SEmin(f) = bestSE;
            continue
        end
        
        fprintf('Fitting: %d %s\n', subjids(f), model);
        
        %get the parameter initial values and limits based on the model to be fit.
        [init_params, lower_limits, upper_limits] = getParamInitialization(model);
        
        if strcmp(model, 'noemo_scram') || strcmp(model, 'noemosticky_scram')
            emoSubset=find(strcmp(emoNames, 'scram'));
        else
            emoSubset=-1;
        end
        
        if generative == 0
            
            %use multiple starting values?
            if multstart == 1
                num_start_pts = 5; % number of initial starting points
                DiffFmOptimal(f,:) = zeros(num_start_pts,1);
                
                opts = optimset('fmincon');
                opts.LargeScale = 'off';
                opts.Algorithm = 'active-set';
                opts.Display = 'none';
                
                %core fitting function -- returns results of length num_start_pts (5) above.
                %These contain fit estimates for each starting point.
                %Then identify the best-fitting output for use in analyses.
                [params, SE, exitflag, xstart] = rmsearch(@(params) TC_minSE(params, subjdata{f}, model, emoSubset), 'fmincon', init_params, ...
                    lower_limits, upper_limits, 'initialsample', num_start_pts, 'options', opts);
                
                SEmin(f)= min(SE);
                
                DiffFmOptimal(f,:) = SE - SEmin(f); % how different are the SSE values for each starting pt from optimal one
                
                %re-run TC_alg for all blocks with optimal parameters
                [totalSqErr, ret_all] = TC_minSE(params(find(SE == min(SE), 1 ),:), subjdata{f}, model, emoSubset);
                
                %[SEbest, PE, exp, std_f, std_s, mn_f, mn_s, Go, NoGo ] = SavePredsFmBest(params(find(SE == min(SE), 1 ),:), subjdata{f}); % save predictions from best run of rmsearch
                
            else
                % use below line if just want to run one starting point (faster and
                % usually not far off from optimal.)
                [params, SE(f), exitflag] = fmincon(@(params) TC_minSE(params, subjdata{f}, model, emoSubset), init_params, [], [], [], [], lower_limits, upper_limits, [], options) ;
                [totalSqErr, ret_all] = TC_minSE(params(find(SE == min(SE), 1 ),:), subjdata{f}, model, emoSubset); %run once with best params... (can we get all this back from fmincon?)
                %SavePredsFmBest(params, subjdata{f});
                
                SEmin(f) = SE;
            end
        else % generative model, used for generating agent-based behavior
            
            SE = TC_minSE(init_params, subjdata{f}, model, emoSubset);
            %RTGene_preds;
            bestFit_all(f, :) = [this_subj s sqrt(SE)] % note these are sqrt of sum!
            
            SEmin(f) = SE;
        end
        
        if multstart==1
            %1 in this vector represents a dummy code for session (just have one session at the moment)
            bestFit_all(f, :) = [subjids(f) 1 params(find(SE == min(SE), 1 ),:) sqrt(SEmin(f))]; % note errors are sqrt of sum!
        else
            bestFit_all(f, :) = [subjids(f) 1 params' sqrt(SE(f))];
        end
        
        
        %for now, create vectors that correspond to what was present before (for checking against vetted code)
        subject = subjids(f);
        PE = reshape([ret_all.rpe], [], 1);
        exp = reshape([ret_all.explore], [], 1);
        std_f = reshape([ret_all.sdShort], [], 1);
        std_s = reshape([ret_all.sdLong], [], 1);
        mn_f = reshape([ret_all.meanShort], [], 1);
        mn_s = reshape([ret_all.meanLong], [], 1);
        Go = reshape([ret_all.go], [], 1);
        NoGo = reshape([ret_all.noGo], [], 1);
        bestFit = bestFit_all(f, :);
        bestSE = SEmin(f);
        save(strcat('../outputs/parameter_mat/modelVars_', num2str(subjids(f)), model), 'subject',  'PE', 'exp', 'std_f', 'std_s', 'mn_f', 'mn_s', 'Go', 'NoGo', 'ret_all', 'bestFit', 'bestSE');
        
        plotSubjData(subjids(f), ret_all, model);
    end
    
    fname_bestFit=sprintf('SubjsSummary_%s.txt', model);
    
    if strcmp(model, 'noemo') || strcmp(model, 'noemo_scram')
        hdr = {'Subject','Session','lambda','explore','alphaG','alphaN','K','nu','ignore','rho','SSE'};
    elseif strcmp(model, 'noemosticky') || strcmp(model, 'noemosticky_scram')
        hdr = {'Subject','Session','lambda','explore','alphaG','alphaN','K','sticky_decay','ignore','rho','SSE'};
    elseif strcmp(model, 'emoexplore')
        hdr = {'Subject','Session','lambda','explore_scram', 'explore_fear', 'explore_happy', 'alphaG','alphaN','K','nu','ignore','rho','SSE'};
    elseif strcmp(model, 'emoexploresticky')
        hdr = {'Subject','Session','lambda','explore_scram', 'explore_fear', 'explore_happy', 'alphaG','alphaN','K','sticky_decay','ignore','rho','SSE'};
    elseif strcmp(model, 'emogonogo')
        hdr = {'Subject','Session','lambda','explore','alphaG_scram','alphaG_fear','alphaG_happy','alphaN_scram','alphaN_fear', 'alphaN_happy', 'K','nu','ignore','rho','SSE'};
    elseif strcmp(model, 'emogonogosticky')
        hdr = {'Subject','Session','lambda','explore','alphaG_scram','alphaG_fear','alphaG_happy','alphaN_scram','alphaN_fear', 'alphaN_happy', 'K','sticky_decay','ignore','rho','SSE'};
    elseif strcmp(model, 'emonogo')
        hdr = {'Subject','Session','lambda','explore','alphaG', 'alphaN_scram', 'alphaN_fear', 'alphaN_happy', 'K','nu','ignore','rho','SSE'};
    elseif strcmp(model, 'emonogosticky')
        hdr = {'Subject','Session','lambda','explore','alphaG', 'alphaN_scram', 'alphaN_fear', 'alphaN_happy', 'K','sticky_decay','ignore','rho','SSE'};
    end

    
    txt=sprintf('%s\t',hdr{:});
    txt(end)='';
    dlmwrite(fname_bestFit,txt,'');
    dlmwrite(fname_bestFit, bestFit_all,'-append','delimiter','\t','precision', '%6.5f');
    
    fname_trn = sprintf('GroupStats_%s.doc', model);
    
    fid_Trn =fopen(fname_trn,'w');
    rSE_Trn_mean = mean(sqrt(SEmin))
    rSE_Trn_std = std(sqrt(SEmin))
    
    fprintf(fid_Trn,'%s \t', 'rSE_Trn_mean = ');
    fprintf(fid_Trn,'%f \n', rSE_Trn_mean);
    fprintf(fid_Trn,'%s \t', 'rSE_Trn_std = ');
    fprintf(fid_Trn,'%f \n', rSE_Trn_std);
    
    fclose(fid_Trn);
    
end

%clear
%matlabpool close