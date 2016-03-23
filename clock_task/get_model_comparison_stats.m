% quantify the advantage of the fixed decay model over others for 
% further analyses, particularly non-model-based
all_fixed_decays = zeros(size(L));
for ct=1:9
    all_fixed_decays(ct,:) = L(3,:);
end
fixed_decay_advantage = all_fixed_decays - L;
fixed_decay_advantage_log_ratio = log(L./all_fixed_decays);

%% check distributions
figure(33); hist(fixed_decay_advantage');
figure(34); hist(fixed_decay_advantage_log_ratio',40);

%% save results
save fixed_decay_advantage_log_ratio fixed_decay_advantage_log_ratio


