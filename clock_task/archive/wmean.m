function R=wmean(X,W,dim)
%weighted mean of X along the dimension dim
%% PERMUTE METHOD
perm_order=[dim 1:dim-1 dim+1:ndims(X)];
X=permute(X, perm_order );
sum1=0;
for i=1:size(X,1)
 sum1 = sum1 + W(i)*X(i,:,:,:,:,:,:,:);
end;
R=sum1/sum(W);
R=ipermute(R,[perm_order]);