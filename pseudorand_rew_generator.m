function b=pseudorand_rew_generator(samples)
% rts = rand(1,trials)*5000;

% %cond = 78;
rts = 1:samples;
rts_real = rts.*(5000./samples);
b.iev_rew = zeros(1,length(rts));
b.dev_rew = zeros(1,length(rts));
for cond = {'DEV' 'IEV', 'QUADUP'}
    for i = 1:length(rts)
        if strcmpi(cond, 'IEV')
            b.iev_rew(i) = RewFunction(rts_real(i),cond);
        elseif strcmpi(cond, 'DEV')
            b.dev_rew(i) = RewFunction(rts_real(i),cond);
        elseif strcmpi(cond, 'QUADUP')
            b.quadup_rew(i) = RewFunction(rts_real(i),cond);
        end
    end
end

b.rts = rts;
b.rts_real = rts_real;
save b;
end