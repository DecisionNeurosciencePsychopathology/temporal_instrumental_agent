function H = calc_shannon_H(p)
    H = sum(-(p(p>0).*(log2(p(p>0)))));
end