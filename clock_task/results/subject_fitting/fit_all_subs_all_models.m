function fit_all_subs_all_models(modelname,sub,subject_param_data,optimal_param_data)

if nargin<3
    num_plots = 1;
    data = subject_param_data;
else
    num_plots = 2;
    data = [subject_param_data; optimal_param_data];
end
factor = [.4 .0003]; %two monitor set up, splits up the figures

for i = 1:num_plots
    %It was getting annoying moving the figs constantly
    set(0,'DefaultFigureUnits','normalized', ...
        'DefaultFigurePosition', [1-i*factor(i),.3,.35,.5]);
    figure(i); clf;
    ret = data(i,:);
    unrew_rts = NaN(size(ret.rt_obs));
    unrew_rts(ret.rew_obs==0) = ret.rt_obs(ret.rew_obs==0);
    disp(['Subject: ',num2str(sub),' ', modelname, ' Params: ',num2str(ret.params);]);
    if i==1
        str = 'Subject Fit';
    else
        str = 'Optimal Fit';
    end
    
    
    
    subplot(3,1,1);
    plot(1:length(ret.rt_obs), ret.rt_obs, 'r');
    hold on;
    plot(1:length(ret.rt_chosen), ret.rt_chosen, 'b');
    hold off;
    title([modelname,' Red: actual RT, Blue: predicted RT ' str]);
    ax2=subplot(3,1,2);
    contourf(1:ret.ntrials, 1:ret.ntimesteps, ret.v_it(1:ret.ntrials,:)'); colorbar('southoutside');hold on;
    scatter(1:ret.ntrials, ret.rt_obs,ret.rew_obs+10, 'r','Filled');
    scatter(1:ret.ntrials, unrew_rts,'b', 'Filled'); hold off;
    title('Value map; red: rewards, blue: omissions');
    colormap(ax2,summer);
    
    if strfind(modelname, 'fixed')
        
        
    elseif  strfind(modelname, 'kalman_processnoise')
        ax3 = subplot(4,1,3);
        contourf(1:ret.ntrials, 1:ret.nbasis, ret.Q_ij(1:ret.ntrials,:)');colorbar('southoutside'); hold on;
         scatter(1:ret.ntrials, ret.rt_obs.*24./500,ret.rew_obs.*24./500+10, 'r','Filled');
    scatter(1:ret.ntrials, unrew_rts.*24./500,'b', 'Filled'); hold off;
         title('Process Noise map; red: rewards,   blue: ommissions');
        colormap(ax3,summer);
         ax3 = subplot(4,1,4);
        contourf(1:ret.ntrials, 1:ret.nbasis, ret.delta_ij(1:ret.ntrials,:)');colorbar('southoutside'); hold on;
         scatter(1:ret.ntrials, ret.rt_obs.*24./500,ret.rew_obs.*24./500+10, 'r','Filled');
    scatter(1:ret.ntrials, unrew_rts.*24./500,'b', 'Filled'); hold off;
         title('Prediction error map; red: rewards,   blue: ommissions');
        colormap(ax3,summer);

    elseif strfind(modelname, 'uv_sum')
        ax3 = subplot(3,1,3);
        contourf(1:ret.ntrials, 1:ret.ntimesteps, ret.u_it(1:ret.ntrials,:)');colorbar('southoutside'); hold on;
        scatter(1:ret.ntrials, ret.rt_obs,ret.rew_obs+10, 'r', 'Filled');
        scatter(1:ret.ntrials, unrew_rts,'b', 'Filled'); hold off;
        title('Uncertainty map; red: rewards,   blue: ommissions');
        colormap(ax3,summer);
    
    elseif strfind(modelname, 'kalman_softmax')
        
    elseif strcmpi(modelname,'qlearning')
        
    end
    
end


waitforbuttonpress;
end