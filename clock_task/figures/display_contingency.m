function display_contingency(cond)
%Choose which contingency you'd like to view
%Options: IEV, DEV, QUADUP
%EX: display_contingency('IEV')

tvec = 1:500; %timesteps

%objective expected value for this function
ev=[];
for val = 1:length(tvec)
        [~,ev(val)] = RewFunction(tvec(val).*10, cond);
end

figure(6); plot(tvec, ev);
title(['Expected value of contingency ' cond]);
ylabel('Value')
xlabel('Time')