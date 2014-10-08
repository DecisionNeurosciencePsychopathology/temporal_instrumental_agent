
% IEV

figure(2); mesh(iev_optimal_params.value_hist)
title('Learning in IEV')
xlabel('Time')
ylabel('Trial')
zlabel('Value')


% DEV

figure(3); mesh(dev_optimal_params.value_hist)
title('Learning in DEV')
zlabel('Value')
ylabel('Trial')
xlabel('Time')