%rescale the truncated basis so that the typical RBF (not truncated) has range 0..1. The current
%implementation is leading to tiny eligibility traces once you multiply by learning rate.
%technically, this needs to be computed for a given basis setup since the AUC=1.0 scaling will depend on the
%number of timesteps. This correction still maintains greater weighting at the edge (i.e., it is proportionate
%to the AUC=1.0) while giving the reasonable behavior of eligibility representing a "proportion of information
%transfer possible" interpretation.

%find the center that is closest to the midpoint of the time interval
%use the value of the regular basis (0..1 scaling) divided by the truncated basis at the max as the adjustment factor
%note that maxauc above will be equal to this correction factor for all RBFs that are not truncated! Thus, we
%are "undoing" the rescaling of these RBFs, while maintaining the correction at the edge.
%[~,indmid]=min(abs(c - median(tvec)));
%trunc_adjust=max(gaussmat(indmid,:))/max(gaussmat_trunc(indmid,:));
%gaussmat_trunc=gaussmat_trunc*trunc_adjust;


%NB: The working combination is: Regular Gaussian RBF for function evaluation, but weight update proceed by truncated
%Gaussian basis and truncated Gaussian spread function. Conceptually, this is like squeezing in all of the variance of
%the update inside the finite interval but then evaluating the function on the regular RBF... to be developed further
%since this remains a bit of a mystery. :)

%note: with truncated gaussian basis (which functions properly), the idea of extending the basis representation outside
%of the finite interval of interest does not apply (i.e., the functions are not defined outside of the bounds).
%gauss mat for all possible timesteps
%lowest_t = tmin - 4*sig; %3 SDs below lowest center.
%highest_t = tmax + 4*sig;

%t_all = round(lowest_t):round(highest_t);

%gaussmat_all = zeros(nbasis,length(t_all));

%for j = 1:nbasis
%    gaussmat_all(j,:) = gaussmf(t_all,[sig c(j)]);
%end
%plot(t_all, gaussmat_all);


%use this to rescale elibility below to maintain constant AUC equivalent to a standard Gaussian membership function
%this leads to eligibility 0-1.0 for eligibility functions within the interval, and > 1.0 max for truncated functions.
%in testing (weightfit.m), this gives the best sampling behavior


%10/16/2014 NB: Even though some matrices such as v_it have columns for discrete timesteps, the
%agent now learns entirely on a continuous time basis by maximizing the value and uncertainty curves over the
%radial basis set and estimating total uncertainty as the definite integral of the uncertainty function over the
%time window of the trial.



% if trial_plots
%     figure(2); clf;
%     %         mesh(1:ntimesteps, 1:ntrials, p_choice);
%     %         xlabel('Response time'); ylabel('Trial'); zlabel('Response probability');hold on;
%     %         scatter3(rt_obs,1:ntrials,  p_chosen,'r*')
% end


%     %choice rule
%
%     % find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)
%
%     % u -- total amount of uncertainty on this trial (starts at 0 and decreases)
%     % u = mean(u_func);
%
%     %use integration to get area under curve of uncertainty
%     %otherwise, our estimate of u is discretized, affecting cost function
%
%     total_u = integral(@(x)rbfeval(x, sigma_ij(i+1,:), c, ones(1,nbasis).*sig), min(tvec), max(tvec));
%     u = total_u/max(tvec); %make scaling similar to original sum? (come back)...
%
%     if u == 0
%         rt_explore = rt_obs(1); %%FEED RTOBS(1) on the first trial so that the fit is not penalized by misfit on first choice
%     else
%         %rt_explore = fminbnd(@(x) -rbfeval(x, sigma_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
%         %leave out any gaussian noise from RT explore because we can't expect subject's data to fit with
%         %random additive noise here.
%         rt_explore = find(u_func==max(u_func), 1); + round(sig_expnoise*randn(1,1)); %return position of first max and add gaussian noise
%
%         if rt_explore < minrt
%             rt_explore = minrt;  %do not allow choice below earliest possible response
%         end
%         if rt_explore > ntimesteps
%             rt_explore = ntimesteps;
%         end
%     end
%
%     %fprintf('trial: %d rt_exploit: %.2f rt_explore: %.2f\n', i, rt_exploit, rt_explore);
%
%     discrim = 0.1;
%
%     %p_explore_i(i) = 1/(1+exp(-discrim.*(u - u_threshold))); %Rasch model with epsilon as difficulty (location) parameter
%
%     %soft classification (explore in proportion to uncertainty)
%     rng(explore_rng_state); %draw from explore/exploit rng
%     choice_rand=rand;
%     explore_rng_state=rng; %save state after random draw above
%
%     %% AD: the GRW exploration leads agent to repeat subject's RT
%     %   solution: skip it for now
%
%     %determine whether to strategic explore versus GRW
%     rng(exptype_rng_seed);
% %     explore_type_rand=rand;
%     exptype_rng_seed=rng;
%
%
%
%
%     %% AD: NB grw exploration is commented out here to avoid circularity with subject's choices
%
%     %rng('shuffle');
%     if i < ntrials
%         if choice_rand < p_explore_i(i)
% %             if (explore_type_rand > k)
%                 %strategic
%                 exptxt='strategic explore';
%                 rt_pred_i(i+1) = rt_explore;
% %             else
% %                 %grw
% %                 exptxt='grw explore';
% %                 rt_grw = rt_obs(i) + grw_step; %use GRW around prior *observed* RT
% %
% %                 %N.B.: Need to have more reasonable GRW near the edge such that it doesn't just oversample min/max
% %                 %e.g., perhaps reflect the GRW if rt(t-1) was already very close to edge and GRW samples in that
% %                 %direction again.
% %                 if rt_grw > max(tvec), rt_grw = max(tvec);
% %                 elseif rt_grw < min(tvec), rt_grw = min(tvec); end
% %                 rt_pred_i(i+1) = rt_grw;
% %             end
%             explore_trials = [explore_trials i+1];
%         else
%             exptxt='exploit';%for graph annotation
%             rt_pred_i(i+1) = rt_exploit;
%             exploit_trials = [exploit_trials i+1]; %predict next RT on the basis of explore/exploit
%         end
%
%         %playing with basis update at the edge
%         %rt_pred_i(i+1) = randi([400,500],1); %force to late times
%
%     end
%%
%
%     %don't compute distance on first trial
%     if (i > 1)
%         %compute deviation between next observed RT and model-predicted RT
%         %use the 2-D Euclidean distance between RT_obs and RT_pred in U and Q space (roughly related to the idea
%         %of bivariate kernel density).
%
%         %normalize V and U vector to unit lengths so that costs do not scale with changes in the absolute
%         %magnitude of these functions... I believe this is right: discuss with AD.
%
%         vnorm=v_it(i,:)/sqrt(sum(v_it(i,:))^2);
%         unorm=u_it(i,:)/sqrt(sum(u_it(i,:))^2);
%
%         %I believe we should use i for this, not i+1 so that we are always comparing the current RT with the
%         %predicted current RT. Technically, the predicted RT was computed on the prior iteration of the loop, but
%         %in trial space, both are current/i.
%         %place RTs into U and Q space
%         d_i(i) = sqrt( (v_it(i, rt_obs(i)) - v_it(i, rt_pred_i(i))).^2 + (u_it(i, rt_obs(i)) - u_it(i, rt_pred_i(i))).^2);
%
%         %d_i(i) = sqrt( (vnorm(rt_obs(i)) - vnorm(rt_pred_i(i))).^2 + (unorm(rt_obs(i)) - unorm(rt_pred_i(i))).^2)
%
%         %weight cost by model-predicted probability of exploration. Use absolute deviation from 0.5 as measure of
%         %uncertainty about exploration/exploitation, and use the deviation as the power (0..1) to which the
%         %distance is exponentiated.
%
%         %.1 power is a hack for now to make the falloff in cost slower near the edges (not just inverted triangle)
%         precision_i(i) = (abs(p_explore_i(i) - 0.5)/0.5)^.1; %absolute symmetric deviation from complete uncertainty
%
%         %r=100;
%         %probs=0:0.01:1;
%         %precision=arrayfun(@(p) abs(p - 0.5)/0.5, probs).^.1;
%         %dr=r.^precision;
%         %plot(probs, dr)
%
%         %note that this multiplying d_i by precision_i gives an inverted triangle where cost = 0 at 0.5
%         %probably want something steeper so that only costs ~0.4-0.6 are severely downweighted
%         d_i(i) = d_i(i)^precision_i(i); %downweight costs where model is ambivalent about explore/exploit
%
%     end
%
%     extra_plots = 0;
%     figure(11);
%     if(extra_plots == 1 && ~all(v_it(i,:)) == 0)
%         gkde2(vertcat(v_it(i,:), u_it(i,:)));
%         %gkde2(vertcat(vnorm, unorm))
%     end
%
%     %compute the expected returns for the observed and predicted choices
%     if (nargin >= 8)
%         [~, ev_obs_i(i)] = RewFunction(rt_obs(i).*10, cond); %multiply by 10 because underlying functions range 0-5000ms
%         [~, ev_pred_i(i)] = RewFunction(rt_obs(i).*10, cond); %multiply by 10 because underlying functions range 0-5000ms
%     end
%
%     verbose=0;
%     if verbose == 1
%        fprintf('Trial: %d, Rew(i): %.2f, Rt(i): %.2f\n', i, rew_obs(i), rt_obs(i));
%        fprintf('w_i,k:    '); fprintf('%.2f ', mu_ij(i,:)); fprintf('\n');
%        fprintf('delta_ij:   '); fprintf('%.2f ', delta_ij(i,:)); fprintf('\n');
%        fprintf('w_i+1,k:  '); fprintf('%.2f ', mu_ij(i+1,:)); fprintf('\n');
%        fprintf('\n');
%
%     end

if trial_plots == 1
    %         figure(1); clf;
    %         subplot(5,2,1)
    %         %plot(tvec,v_func);
    %         scatter(rt_pred_i(1:i),rew_obs(1:i)); axis([1 500 0 350]);
    %         hold on;
    %         plot(rt_pred_i(i),rew_obs(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
    %         hold off;
    %         subplot(5,2,2)
    %         plot(tvec,v_func); xlim([-1 ntimesteps+1]);
    %         ylabel('value')
    %         subplot(5,2,3)
    %
    % %         bar(c, mu_ij(i,:));
    % %         ylabel('basis function heights');
    %         plot(tvec,v_jt);
    %         ylabel('temporal basis function')
    % %         title(sprintf('trial # = %i', h)); %
    %                 xlabel('time(ms)')
    %                 ylabel('reward value')
    %
    %         subplot(5,2,4)
    %         plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
    %         xlabel('time (centiseconds)')
    %         ylabel('uncertainty')
    %
    %         subplot(5,2,5)
    %         barh(sigmoid);
    %         xlabel('p(explore)'); axis([-.1 1.1 0 2]);
    %         subplot(5,2,6)
    %         %barh(alpha), axis([0 .2 0 2]);
    %         %xlabel('learning rate')
    %         subplot(5,2,7)
    %         barh(sig_spread), axis([0.0 0.5 0 2]);
    %         xlabel('decay')
    %         subplot(5,2,8)
    %         barh(epsilon), axis([-0.5 0 0 2]);
    %         xlabel('strategic exploration')
    %
    %         subplot(5,2,9)
    %         barh(u) %, axis([0 1000 0 2]);
    %         xlabel('mean uncertainty')
    %         %         pause(0.1);
    %         subplot(5,2,10)
    %         plot(1:ntrials,rt_pred_i, 'k');
    %         ylabel('rt by trial'); axis([1 ntrials -5 505]);
    %
    %         figure(1); clf;
    %         set(gca,'FontSize',18);
    %         subplot(3,2,1);
    %         title('Choice history');
    %         %plot(tvec,v_func);
    %         scatter(rt_obs(1:i),rew_obs(1:i)); axis([1 ntimesteps 0 350]);
    %         hold on;
    %         plot(rt_obs(i),rew_obs(i),'r*','MarkerSize',20);  axis([1 ntimesteps 0 350]);
    %         hold off;
    %         subplot(3,2,2)
    %         title('Learned value');
    %         plot(tvec,v_func); xlim([-1 ntimesteps+1]);
    %         ylabel('expected value')
    % %         bar(c, mu_ij(i,:));
    % %         ylabel('basis function heights');
    %         %title('basis function values');
    %         %plot(tvec,v_jt);
    %         %ylabel('temporal basis function')
    % %         title(sprintf('trial # = %i', h)); %
    %
    %         subplot(3,2,3);
    %         scatter(rt_pred_i(1:i),rew_obs(1:i)); axis([1 ntimesteps 0 350]);
    %         %         text(20, max(rew_obs), exptxt);
    %         hold on;
    %         plot(rt_pred_i(i),rew_obs(i),'r*','MarkerSize',20);  axis([1 ntimesteps 0 350]);
    %         hold off;
    %
    %         subplot(3,2,4);
    %         plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
    %         xlabel('time (centiseconds)')
    %         ylabel('uncertainty')
    %
    %         subplot(3,2,5);
    %         %eligibility trace
    %         title('eligibility trace');
    %         %elig_plot = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 1);
    %         %plot(tvec, elig_plot);
    %         plot(tvec, elig);
    %         xlabel('time(centiseconds)')
    %         ylabel('eligibility')
    %
    %
    %         subplot(3,2,6);
    %         plot(1:length(rt_obs), rt_obs, 'r');
    %         hold on;
    %         plot(1:length(rts_pred_exploit), rts_pred_exploit, 'b');
    % %         plot(1:length(rts_pred_explore), rts_pred_explore, 'k');
    % %         plot(1:length(rts_pred_explore), rts_uv_pred, 'g');
    %         hold off;
    %
    %         %figure(2); clf;
    %         %plot(tvec, u_func);
    %         %hold on;
    %         %plot(c, e_ij(i,:))
    %         %plot(c, e_ij(1:i,:)')
    %         %bar(c, sigma_ij(i,:))
    %
    %         drawnow update;
    %         mov(i) = getframe(gcf);
    %     end
    %     disp([i rt_pred_i(i) rew_obs(i) sum(v_func)])
end

%ret.cost_explore = -sum((rt_pred_i(explore_trials) - rt_obs(explore_trials)).^2);
%ret.cost_exploit = -sum((rt_pred_i(exploit_trials) - rt_obs(exploit_trials)).^2);


%%OLD APPROACH TO CHOOSING RT EXPLOIT WITH CONTINUOUS BASIS

%DANGER: fminbnd is giving local minimum solution that clearly
%misses the max value. Could switch to something more comprehensive
%like rmsearch, but for now revert to discrete approximation
%rt_exploit = fminbnd(@(x) -rbfeval(x, mu_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
%figure(2);
%vfunc = rbfeval(0:500, mu_ij(i+1,:), c, ones(1,nbasis).*sig);
%plot(0:500, vfunc);

% find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)


%OLD APPROACH TO COMPUTING AUC OF UNCERTAINTY
% u -- total amount of uncertainty on this trial (starts at 0 and decreases)
% u = mean(u_func);

%use integration to get area under curve of uncertainty
%otherwise, our estimate of u is discretized, affecting cost function
%total_u = integral(@(x)rbfeval(x, sigma_ij(i+1,:), c, ones(1,nbasis).*sig), min(tvec), max(tvec));
%u = total_u/max(tvec); %make scaling closer to value per timestep