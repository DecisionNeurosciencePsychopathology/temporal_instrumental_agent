function [cost]=correctrbf(weights, tvec, centers, widths, trunc)
    funcapprox=rbfeval(tvec, weights, centers, widths, trunc);
    cost=sum((funcapprox - 1).^2); %variation from unit magnitude
end
