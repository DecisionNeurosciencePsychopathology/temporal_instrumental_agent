%Determine who wins the most via posterior from groupBMC
[M,I] =max(posterior.r);

for i = 1:length(modelnames)
    m(1,i) = sum(I==i);
    subjs{1,i} = find(I==i);
end