function checkDecayDistrib
file_path = '/Volumes/bek/vba_results/fixed_decay/';
s = getPosteriorVariable(file_path,'muTheta');
gamma = (1./(1+exp(-s.muTheta(:,2))))./1;%Transform this!
getParamterDistribtion(gamma);
corr(gamma,s.outF)




function getParamterDistribtion(var)
    histogram(var)



function var_struct = getPosteriorVariable(path_def,vars_of_interest)
    file_pointers = glob([path_def, '*.mat']);
        
    parfor i = 1:length(file_pointers)
        data = load(file_pointers{i});
        muTheta(i,:) = data.posterior.(vars_of_interest);
        outF(i,:) = data.out.F;
    end
    var_struct.muTheta = muTheta;
    var_struct.outF = outF;
