function res=rbfeval(val, weights, centers, widths)

%is there a way to vectorize?
res=nan(1,length(val));
for v=1:length(val)
    cdiff = centers - val(v);
    res(v) = sum(weights .* exp(-cdiff.^2./(2.*widths.^2)));
end

end