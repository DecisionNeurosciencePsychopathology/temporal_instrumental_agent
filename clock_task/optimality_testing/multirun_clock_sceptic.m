function [totcost, costs] = multirun_clock_sceptic(params, agent, runarray)
%use multiple runs of data to identify optimal parameters for sceptic agent
nruns = length(runarray);
costs=NaN(nruns, 1);

reversal = 0;
for i = 1:nruns
    %args{1,9}.lookup = args{1,9}.lookup(:,vperm(i,:)); %Permute
    thiscall={params, agent, runarray(i).seeds, runarray(i), agent.ntrials, agent.nbasis, agent.ntimesteps, reversal};
    costs(i) = clock_sceptic_agent(thiscall{:});
end

totcost=sum(costs);
%fprintf('totcost: %.2f, pars: %s\n', totcost, num2str(params));
end

