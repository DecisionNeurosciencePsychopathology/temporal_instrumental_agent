function [cost, yhat, centers]=rbffit(params, x, y, model, sds, margin_offset)
    if nargin < 4, model = 1; end    
    if nargin < 5, sds = 0.3; end
    if nargin < 6, margin_offset = 0; end
    
    %default to evenly spaced centers
    xmin=min(x) - margin_offset; xmax=max(x) + margin_offset;
    
    if model == 1
        %weights free, sd fixed, centers fixed
        weights=params;
        nbasis=length(weights);
        centers=xmin:(xmax-xmin)/(nbasis-1):xmax;
    elseif model == 2
        %weights free, sd free, centers fixed
        sds=params(1);
        weights=params(2:end);
        nbasis=length(weights);
        centers=xmin:(xmax-xmin)/(nbasis-1):xmax;
    elseif model == 3
        %weights free, sd free, centers free
        sds=params(1);
        centers=params(2:((length(params)-1)/2)+1);
        weights=params(((length(params)-1)/2)+2:end);
    end
    
    %fprintf('weights: %s\n', num2str(weights));
    %fprintf('centers: %s\n', num2str(centers));
    %fprintf('sd: %.3f\n', sds);

    %evaluate the current RBF at all x values to get predicted y values
    yhat = rbfeval(x, weights, centers, sds)';
    cost = sum((y - yhat).^2);
    %fprintf('param values: %s\n', num2str(params));
    fprintf('cost: %.3f\n', cost);
end
