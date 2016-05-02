function sceptic_movie(costs, a, ret, pars, optmat, outname)
writerObj = VideoWriter(outname); % Name it.
writerObj.FrameRate = 8; % How many frames per second.open(writerObj);
open(writerObj); %Need to open writeerObj first

figure('units','pixels','position',[0 0 1000 700]);

for rep = 1:length(costs)
    ntrials = length(ret{rep}.rts);
    t = 1:ret{rep}.ntimesteps;
    %pars_rep = pars(rep,:);
    
    %inv logit gamma
    if strcmpi(a.name, 'fixedLR_decay'),
        pars(rep,strcmpi(a.parnames,'gamma')) = 1/(1+exp(-pars(rep,strcmpi(a.parnames,'gamma')))); %exponentiate to convert to 0-1 scaling
    end
    
    bestcost=max(optmat(rep).ev)*ntrials;
    %clunky!
    parnames=sprintf('Rep: %d, Earned: %.1f, Best: %.1f, ', rep, -1*costs(rep), bestcost);
    for kk = 1:size(pars,2)
        parnames = [parnames, a.parnames{kk}, ' = ', sprintf('%.3f', pars(rep,kk))];
        if kk < size(pars,2), parnames=[parnames, ', ']; end
    end
    
    for i = 2:ntrials
        gcf; th=suptitle(parnames); set(th,'interpreter','none');
        subplot(3,3,1:3);
        title('True expected value (red), Vfinal (green), and choices (blue stars)');
        plot(1:a.ntimesteps, optmat(rep).ev,'r','LineWidth',2);
        hold on;
        plot(ret{rep}.rts(1:i) + rand(1,i) - 0.5,optmat(rep).ev(ret{rep}.rts(1:i)) + 1,'*b');
        plot(1:a.ntimesteps, ret{rep}.vfinal_it(i,:), 'g', 'LineWidth', 2);
        hold off;
        %change in value function from i-1 to i
        subplot(3,3,4:6);
        title('V final i (green) versus i-1 (blue)');
        plot(t, ret{rep}.vfinal_it(i,:),'Color',[0 .5 0],'LineWidth',2);
        hold on;
        plot(t, ret{rep}.vfinal_it(i-1,:),'b','LineWidth',2);
        text(ret{rep}.rts(i), (max(ret{rep}.vfinal_it(i,:)) - min(ret{rep}.vfinal_it(i,:)))*0.5, sprintf('%.1f', ret{rep}.rew_i(i)), 'HorizontalAlignment', 'center');
        hold off;
        subplot(3,3,7);
        title('RT histogram');
        hist(ret{rep}.rts(1:i));
        subplot(3,3,8:9);
        title('p choice');
        plot(t, ret{rep}.p_choice(i,:));
        frame = getframe(gcf);
        
        writeVideo(writerObj, frame);
    end
end

close(writerObj);
close(gcf);
end