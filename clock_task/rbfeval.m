function res=rbfeval(val, weights, centers, widths, trunc)

computegaussmat=1;
if nargin == 3
    %gaussmat is passed in directly as third parameter so does not need to be recomputed
    gaussmat = centers;
    computegaussmat = 0;
elseif nargin < 5
    trunc=[-Inf Inf];
end

pdfs=[];

if ~computegaussmat
    %compute gaussmat for each basis j and timestep t by multiplying weights by basis
    g_jt=weights'*ones(1,size(gaussmat,2)) .* gaussmat; %use vector outer product to replicate weight vector
    
    res=sum(g_jt);
elseif trunc(1) ~= -Inf
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
    %typical non-truncated basis.
    %better to vectorize for speed as above (this is taking up a ton of time according to the profiler)
    %this version is slow! (kept here for other algs that haven't implemented the update)
    res=nan(1,length(val));
    
    for v=1:length(val)
        cdiff = centers - val(v);
        res(v) = sum(weights .* exp(-cdiff.^2./(2.*widths.^2)));
    end
    
end

end