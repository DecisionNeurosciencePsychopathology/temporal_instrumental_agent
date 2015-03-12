function res=rbfeval(val, weights, centers, widths, trunc)

if nargin < 5
    trunc=[-Inf Inf];
end

pdfs=[];
if trunc(1) ~= -Inf
    %switch to a truncated gaussian basis so that pdf respects probabilities within 0-500 bounds
    for i=1:length(centers)
        ntemp=makedist('Normal','mu',centers(i),'sigma',widths(i));
        pdfs{i} = truncate(ntemp, trunc(1), trunc(2));
    end
    
    pdfeval = cellfun(@(p) pdf(p, val), pdfs, 'UniformOutput', false);
    res=zeros(1,length(val));
    for i=1:length(pdfeval)
        res=res + weights(i) .* pdfeval{i};
    end

else
    %old update
    %is there a way to vectorize?
    res=nan(1,length(val));
    for v=1:length(val)
        cdiff = centers - val(v);
        res(v) = sum(weights .* exp(-cdiff.^2./(2.*widths.^2)));
    end

end

end