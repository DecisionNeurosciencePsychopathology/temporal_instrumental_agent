
% IEV

figure(1); mesh(iev_fitted_higher_range.value_hist)
title('Learning in IEV')
xlabel('Time')
ylabel('Trial')
zlabel('Value')


% DEV

figure(2); mesh(dev_fitted_higher_range.value_hist)
zlabel('Value')
ylabel('Trial')
xlabel('Time')