function [totcost, costs, seeds] = multirun_clock_sceptic(params, agent, runarray)
%use multiple runs of data to identify optimal parameters for logistic
%operator.
rng(agent.runseed);
nruns = length(runarray);
seeds=randi([1 500], nruns, 4);
costs=NaN(nruns, 1);

reversal = 0;
for i = 1:nruns
    %args{1,9}.lookup = args{1,9}.lookup(:,vperm(i,:)); %Permute
    thiscall={params, agent.name, seeds(i,:), runarray(i), agent.ntrials, agent.nbasis, agent.ntimesteps, reversal};
    costs(i) = clock_sceptic_agent(thiscall{:});
end

totcost=sum(costs);
end
