allsubjs = dir('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/*.csv');
nsubjs = length(allsubjs);
ids = NaN(nsubjs, 1);
% for i = 1:length(allsubjs)
%     ids(i) = str2double(regexprep(allsubjs(i).name,'fMRIEmoClock_(\d+)_tc_tcExport.csv','$1'));
%     fprintf('processing id: %d\n', ids(i));
%     %only run if no previous file
%     if ~(exist(sprintf('results/%d_tc_vba_fit.mat', ids(i)), 'file') == 2)
%         [posterior, out] = clock_tc_vba(ids(i), 0); %suppress figure
%     else
%         fprintf('   exists: results/%d_tc_vba_fit.mat\n', ids(i));
%     end
%     
%     if ~(exist(sprintf('results/%d_klambda_vba_fit.mat', ids(i)), 'file') == 2)
%         [posterior, out] = clock_klambda_vba(ids(i), 0); %suppress figure
%     else
%         fprintf('   exists: results/%d_klambda_vba_fit.mat\n', ids(i));
%     end
% end

fullmodel = 'K_Lambda_Nu_AlphaG_AlphaN_Rho_Epsilon';
split = strsplit(fullmodel, '_');
models = {};
for m = 1:length(split)
   models{m} = strjoin(split(1:m), '_'); 
end

clear m

%re-run full model using epsilon unif from 0..10000 instead of 0..80000
models={'K_Lambda_Nu_AlphaG_AlphaN_Rho_Epsilon'};

posteriors = cell(nsubjs, length(models));
outputs = cell(nsubjs, length(models));

parfor i = 1:nsubjs
    ids(i) = str2double(regexprep(allsubjs(i).name,'fMRIEmoClock_(\d+)_tc_tcExport.csv','$1'));
    fprintf('processing id: %d\n', ids(i));
    %only run if no previous file
    
    ipost = cell(1, length(models));
    iout = cell(1, length(models));
    for j = 1:length(models)
        [ipost{j}, iout{j}] = clock_tc_components_vba(ids(i), 0, models{j}); %suppress figure
    end
    posteriors(i,:) = ipost;
    outputs(i,:) = iout;
    
end

%save('allfranktc2.mat', 'posteriors', 'outputs', 'models', 'ids', '-v7.3');

%RFX BMC
% models={'tc', 'klambda'};
% L = NaN(length(models), length(ids));
% for i = 1:length(models)
%     for j = 1:length(ids)
%         load(sprintf('results/%d_%s_vba_fit.mat', ids(j), models{i}));
%         L(i,j) = out.F;
%     end
% end
% 
% load '/Users/michael/Downloads/L_K_lambda.mat';
% 
% Lall = [L; L_K_lambda];
% [posterior,out] = VBA_groupBMC(Lall);
%load('allfranktc.mat')
%new frank incremental variant BMC
logEvidence = NaN(size(outputs));
parameters = NaN([size(outputs), 7]);
for i = 1:nsubjs
    for j = 1:length(models)
        logEvidence(i,j) = outputs{i,j}.F;
        %if ismember('transformed', fields(posteriors{i,j})), parameters(i,j,1:j) = struct2array(posteriors{i,j}.transformed); end;
        if ismember('transformed', fields(posteriors{i,j})), parameters(i,j,1:7) = struct2array(posteriors{i,j}.transformed); end;
    end
end
% 
% nanmean(parameters, 2)
% 
% 
% %bmc expects it to be models x evidence -- transpose
% logEvidence = logEvidence';
% save('tc_logevidence.mat', 'logEvidence', 'ids', 'models');
% 
% [BMCposterior,BMCout] = VBA_groupBMC(logEvidence);