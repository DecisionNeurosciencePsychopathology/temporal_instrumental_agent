allsubjs = dir('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/*.csv');
ids = NaN(length(allsubjs), 1);
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
for i = 1:length(split)
   models{i} = strjoin(split(1:i), '_'); 
end

posteriors = cell(length(allsubjs), length(models));
outputs = cell(length(allsubjs), length(models));

parfor i = 1:length(allsubjs)
    ids(i) = str2double(regexprep(allsubjs(i).name,'fMRIEmoClock_(\d+)_tc_tcExport.csv','$1'));
    fprintf('processing id: %d\n', ids(i));
    %only run if no previous file
    
    ipost = cell(length(models), 1);
    iout = cell(length(models), 1);
    for j = 1:length(models)
        [ipost{j}, iout{j}] = clock_tc_components_vba(ids(i), 0, models{j}); %suppress figure
    end
    posteriors{i,:} = ipost;
    outputs{i,:} = iout;
    
end

save('allfranktc.mat', 'posteriors', 'outputs', 'models');

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