function [mov] = plot_sceptic(ret)
%plot sceptic outputs on the basis of the return structure ret
%rather than call plots within trial loop, this function is designed to replay agent behavior
%and return a mov struct compatible with movie generaiton



%initialize movie storage
mov=repmat(struct('cdata', [], 'colormap', []), ntrials,1);




figure(1); clf;
set(gca,'FontSize',18);
subplot(4,2,1);
title('Choice history');
%plot(tvec,v_func);
scatter(rts(1:i),rew_i(1:i)); axis([1 500 0 350]);
text(20, max(rew_i), exptxt);
hold on;
plot(rts(i),rew_i(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
hold off;
subplot(4,2,2)
title('Learned value');
plot(tvec,v_func); xlim([-1 ntimesteps+1]);
ylabel('expected value')
subplot(4,2,3);

%         %eligibility trace
%         title('eligibility trace');
%         %elig_plot = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 1);
%         %plot(tvec, elig_plot);
%         plot(tvec, elig);
% %         bar(c, mu_ij(i,:));
% %         ylabel('basis function heights');
%         %title('basis function values');
%         %plot(tvec,v_jt);
%         %ylabel('temporal basis function')
% %         title(sprintf('trial # = %i', h)); %
%         xlabel('time(centiseconds)')
%         ylabel('eligibility')


subplot(4,2,3);
title('Decay')
plot(k_top_plot)
ylabel('Decay, despair(red)');
hold on;
plot(despair, 'r');
hold off

subplot(4,2,4);
plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
xlabel('time (centiseconds)')
ylabel('uncertainty')

%figure(2); clf;
%plot(tvec, u_func);
%hold on;
%plot(c, e_ij(i,:))
%plot(c, e_ij(1:i,:)')
%bar(c, sigma_ij(i,:))

subplot(4,2,5);
title('RT history');
plot(1:ntrials, rts(1:ntrials));
xlim([0 ntrials]);
ylabel('RT(ms)');


if uvsum==1
    subplot(4,2,6);
    title('UV');
    plot(tvec,uv); xlim([-1 ntimesteps+1]);
    ylabel('UV');
    
    
    subplot(4,2,7);
    title('Prediction Error');
    plot(delta_func); xlim([1 ntrials]);
    hold on
    plot(v_max_to_plot, 'r');
    ylabel('Prediction error');
    
    subplot(4,2,8);
    title('UV-Softmax');
    plot(uv_softmax); xlim([-1 ntimesteps+1]);
    ylabel('UV-Softmax');
end


drawnow update;
mov(i) = getframe(gcf);

end