%% calculate probability of k consecutive "corrects" in a sequence of n trials

k = 6; % e.g. hexanacci number for k=6
n = 50; % number of trials +2


% 
% f=2;
% for i=3:k
%     f = f+fibonacci(n+2);
% end
% 
% f=double(f);
% p = f/(4^n)


%% run 10000 permutations

n_perm = 1000000;

parfor i = 1:n_perm
    rawseq = rand(1,n).*4
    seq = round(rawseq)
    %% find runs
    a = (diff(seq));
    v = find(conv(double(a==0),ones(1,k),'valid')==k); %// find k zeros
    if isempty(v)
        count(i) = 0;
    else
    v = v([true diff(v)>n]); %// remove similar indices, indicating n+1, n+2... zeros
    count(i) = length(v)
    end
    %% find climbs
    b=diff(rawseq);
    rt_increases = b<0;
    climbs = find(conv(double(rt_increases==0),ones(1,k),'valid')==k); %// find k zeros
    if isempty(climbs)
        climb_count(i) = 0;
    else
    climbs = climbs([true diff(climbs)>n]); %// remove similar indices, indicating n+1, n+2... zeros
    climb_count(i) = length(climbs)
    end
end
p_climb = mean(climb_count)
p = mean(count)
figure(99); clf; hist(count);