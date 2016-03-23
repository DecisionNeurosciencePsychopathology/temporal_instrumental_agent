files = glob('/Volumes/bek/vba_results/sigma_volatility_variants/precision/*.mat');

for i = 1:length(files)
    load(files{i})
    L(i) = out.F;
    gamma(i) = 
end