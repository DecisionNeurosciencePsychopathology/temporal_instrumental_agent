X = subj_data.id11280.rtsIEV';
%X = reshape(X,3,50);
[idx,C] = kmeans(X,2);
figure;
plot(X(idx==1,1),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),'b.','MarkerSize',12)
plot(C(:,1),'kx',...
     'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'