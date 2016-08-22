%% plot BMC

figure(1); clf;
subplot(2,1,1)
barh(out.Ef, 'YTickLabel',modelnames); 
axis([0 1 0 10])

subplot(2,1,2)

barh(out.families.Ef, 'r')
axis([0 1 0 4])

