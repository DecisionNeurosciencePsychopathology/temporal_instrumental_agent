function y=winsorising(varargin)
% Winsorising extreme values.
% Winsorising or Winsorization is the transformation of statistics by
% limiting extreme values in the statistical data to reduce the effect of
% possibly spurious outliers. It is named after the engineer-turned-biostatistician 
% Charles P. Winsor (1895–1951). The effect is the same as clipping in 
% signal processing. The distribution of many statistics can be heavily 
% influenced by outliers. A typical strategy is to set all outliers to a
% specified percentile of the data; for example, a 90% Winsorisation would
% see all data below the 5th percentile set to the 5th percentile, and data
% above the 95th percentile set to the 95th percentile. Winsorised
% estimators are usually more robust to outliers than their more standard
% forms, although there are alternatives, such as trimming, that will
% achieve a similar effect. Note that Winsorizing is not equivalent to
% simply excluding data, which is a simpler procedure, called trimming. 
% In a trimmed estimator, the extreme values are discarded; in a Winsorized
% estimator, the extreme values are instead replaced by certain percentiles
% (the trimmed minimum and maximum).  
% Thus a Winsorized mean is not the same as a truncated mean. For instance,
% the 5% trimmed mean is the average of the 5th to 95th percentile of the
% data, while the 90% Winsorised mean sets the bottom 5% to the 5th
% percentile, the top 5% to the 95th percentile, and then averages the data   
% 
% By itself, WINSORISING runs a demo
%
% Syntax: 	Y=winsorising(X,W)
%           X - data matrix or vector
%           For vectors, WINSORISING(X) is the winsorized X array. 
%           For matrices, WINSORISING(X) is a matrix containing the 
%           winsorized element from each column. 
%
%           W - Amount of winsoritazion (90 by default). If you set W=90 this 
%           means that the remaing 10% (0-5th percentile and 95-100th
%           percentile) will be substituted.
%
%      Example: 
%
% x=[92 19 101 58 103 91 26 78 10 13 0 101 86 85 15 89 89 25 2 41];
%
%           Calling on Matlab the function: winsorizing(x)
%
%           Answer is:
%
% y=92 19 101 58 102 91 26 78 10 13 1 101 86 85 15 89 89 25 2 41
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2011). WINSORISING: WINSORISING Data
% http://www.mathworks.com/matlabcentral/fileexchange/

%Input Error handling
args=cell(varargin);
nu=numel(args);
if nu>2 
    error('stats:winsorising:TooMuchInputs','Max two inputs are required.');
else
    default.values = {[92 19 101 58 103 91 26 78 10 13 0 101 86 85 15 89 89 25 2 41],90}; %default values
    default.values(1:nu) = args;
    [x w] = deal(default.values{:});
    if ~all(isnumeric(x(:))) || ~all(isfinite(x(:)))
        error('The X values must be numeric and finite')
    end
    if nu==2 %if b was given
        if ~isscalar(w) || ~isnumeric(w) || ~isfinite(w)
            error('Warning: W value must be scalar, numeric and finite');
        end
        if w <= 0 || w >= 100 %check if b is between 0 and 100
            error('Warning: W must be comprised between 0 and 100.')
        end        
    end
    rw=(100-w)/2;
    b=[rw 100-rw];
end
clear args default nu rw w

lb=prctile(x,min(b)); ub=prctile(x,max(b)); %set lower and upper bound

if isvector(x)
    y=x; y(y<lb)=lb; y(y>ub)=ub; %winsorising using logical indexing
else
    y=zeros(size(x));
    for I=1:size(x,2)
        c=x(:,I);
        c(c<lb(I))=lb(I); c(c>ub(I))=ub(I);
        y(:,I)=c;
    end
end