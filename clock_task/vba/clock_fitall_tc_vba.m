ncpus=getenv('matlab_cpus');
if strcmpi(ncpus, '')
    ncpus=40;
    fprintf('defaulting to 40 cpus because matlab_cpus not set\n');
else
    ncpus=str2double(ncpus);
end

basedir='/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/subjects';

allsubjs = dir([basedir, '/*.csv']);
nsubjs = length(allsubjs);
ids = NaN(nsubjs, 1);

fullmodel = 'K_Lambda_Nu_AlphaG_AlphaN_Rho_Epsilon';
split = strsplit(fullmodel, '_');
models = cell(1, length(split));
for m = 1:length(split)
    models{m} = strjoin(split(1:m), '_');
end

clear m

%uncomment to run only the full model
%models={'K_Lambda_Nu_AlphaG_AlphaN_Rho_Epsilon'};

posteriors = cell(nsubjs, length(models));
outputs = cell(nsubjs, length(models));

poolobj = parpool('local', ncpus);

try
    parfor i = 1:nsubjs
        ids(i) = str2double(regexprep(allsubjs(i).name,'fMRIEmoClock_(\d+)_tc_tcExport.csv','$1'));
        fprintf('processing id: %d\n', ids(i));
        ipost = cell(1, length(models));
        iout = cell(1, length(models));
        for j = 1:length(models)
            [ipost{j}, iout{j}] = clock_tc_components_vba(allsubjs(i).name, 0, models{j}, basedir); %suppress figure (0)
        end
        posteriors(i,:) = ipost;
        outputs(i,:) = iout;
    end
catch err
    disp('error in optimization. killing parpool');
    delete(poolobj);
    rethrow(err);
end

delete(poolobj);

save('allfranktc.mat', 'posteriors', 'outputs', 'models', 'ids', '-v7.3');

%RFX BMC
%load('allfranktc.mat')
%new frank incremental variant BMC
logEvidence = NaN(size(outputs));
parameters = NaN([size(outputs), 7]);
rawparameters = NaN([size(outputs), 7]);
for i = 1:nsubjs
    for j = 1:length(models)
        logEvidence(i,j) = outputs{i,j}.F;
        rawpars = [posteriors{i,j}.muTheta, posteriors{i,j}.muPhi];
        if ismember('transformed', fields(posteriors{i,j}))
            parvec = struct2array(posteriors{i,j}.transformed);
            parameters(i,j,1:length(parvec)) = parvec;
        end
    end
end

% nanmean(parameters, 2)

% %bmc expects it to be models x evidence/subjects -- transpose
logEvidence = logEvidence';
save('tc_logevidence.mat', 'logEvidence', 'ids', 'models', 'parameters', 'rawparameters');
%
%[BMCposterior,BMCout] = VBA_groupBMC(logEvidence);
