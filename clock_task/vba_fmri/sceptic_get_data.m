function [data, y, u] = sceptic_get_data(data_file, n_steps)
  data = readtable(data_file,'Delimiter',',','ReadVariableNames',true);
  range_RT = 400; %10ms representation
  
  n_t = size(data,1); %number of rows
  trialsToFit = 1:n_t; %fit all by default
  n_runs = n_t/50; %50 trials per run
  
  rtrnd = round(data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
  rtrnd(rtrnd==0)=1;

  %% compute multinomial response
  y = zeros(n_steps, length(trialsToFit));
  for i = 1:length(trialsToFit)
    y(rtrnd(i), i) = 1;
  end
  
  %% Inputs
  u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'];
  u = [zeros(size(u,1),1) u(:,1:end-1)]; %right shift inputs to get t versus t+1 correct in inputs versus predictions
  
end
