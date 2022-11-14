function [data, y, u, so] = sceptic_get_data(data_file, so)
if(strcmpi(so.dataset,'explore'))
    load(data_file,'-mat','order','orderfmt')
    order=order(~cellfun('isempty',order));
    outx=cellfun(@cell2table, order, 'UniformOutput', false);
    data=vertcat(outx{:});
    data.Properties.VariableNames=orderfmt;
    writetable(data,strrep(data_file,'.mat','.csv'))
else
    data = readtable(data_file,'Delimiter',',','ReadVariableNames',true);
end

% determine block number by looking at change in rewFunc and emotion
data{:,'block_number'} = 0;
curval = strcat(data{1,'rewFunc'}, data{1,'emotion'});
curit = 1;
for (i = 1:size(data,1))
  testval = strcat(data{i,'rewFunc'}, data{i,'emotion'});
  if (~strcmpi(curval, testval))
    curval = testval;
    curit = curit + 1;
  end
  data{i,'block_number'} = curit;
end

% determine trials per run from data since it might vary (unlikely, but happens in bsocial)
boundaries = find(ischange(data{:,'block_number'})) - 1;
so.trials_per_run = boundaries(1);

n_t = size(data,1); %number of rows
trials_to_fit = 1:n_t; %fit all trials by default
n_runs = n_t/so.trials_per_run; %usually 50 trials per run

% discretize RTs down to time bins
rtrnd = round(data{trials_to_fit,'rt'} / so.trial_length * so.ntimesteps)';
rtrnd(rtrnd==0) = 1; %if rounding results in a time bin of 0, round up to 1
rtrnd(rtrnd > so.ntimesteps) = so.ntimesteps; %if rounding results in a time bin outside of the valid range (e.g., 41), round to max

% compute multinomial response
y = zeros(so.ntimesteps, length(trials_to_fit));
for i = 1:length(trials_to_fit)
    y(rtrnd(i), i) = 1;
end

% inputs
u = [(data{trials_to_fit, 'rt'} / so.trial_length * so.ntimesteps)'; data{trials_to_fit, 'score'}'];
u = [zeros(size(u,1),1) u(:,1:end-1)]; %right shift inputs to get t versus t+1 correct in inputs versus predictions

%add run boundary marker for u resetting
u = vertcat(u, [0; diff(data.run)]');

end
