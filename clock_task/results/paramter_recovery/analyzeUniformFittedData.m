function analyzeUniformFittedData(org_data,fitted_data)
%Reutrn scatter plots of fitted data from rmsearch and residuals
%ex analyzeUniformFittedData(uniform_test_data,uniform_fitted_results)

model_names = fieldnames(org_data);

for i = 1:length(model_names)
    org_params = org_data.(model_names{i}).params;
    %We need to remove the column of data containing Beta from the org
    %params
    idx=org_params(1,:)~=.1;
    
    %Special qlearning case!
    if strcmpi(model_names{i},'qlearning')
        idx=logical([0 1 1 0]);
    elseif strcmpi(model_names{i},'franktc_fixed')
        idx=logical([0 1 1 1 0 0 1]);
    end
    
    %Select only the fitted params
    org_params = org_params(:,idx);
    fitted_params = fitted_data.(model_names{i}).best_parameters;
    
    n_params = size(fitted_params,2); %Number of params
    n_subjs = length(fitted_params);
    n_plot_panels = n_params*2;
    sub=1; %just a place holder
    [~,p_names]=getOrginalParams(model_names{i});
    
    
    %To make figures look better rmv underscores from param names and model
    %name
    if cell2mat(strfind(p_names, '_'))
        
        p_names = cellfun(@strrep, p_names, repmat({'_'},length(p_names),1), repmat({' '},length(p_names),1), 'UniformOutput', false);
    end
    
    if strfind(model_names{i}, '_')
        model_names{i} = strrep(model_names{i},'_',' ' );
    end
    
    
    
    
    %Create a figure for every parameter
    h=figure(i);
    clf;
    %title(model_names{i}) %Need to fix title?
    ct=0;
    x_text = 0; %x coordinate of r2 value on figure
    
    %Set up color map
%     b = [ 0.0 1.0 1.0;
% 1.0 0.5 0.0;
% 1.0 0.0 0.0];
    colormap default;
    for j = 1:n_params
        if j>1,ct=j-1;end
        subplot(n_params,2,j+ct)
        x=org_params(:,j);
        y=fitted_params(:,j);
        y_text = max(y); %y coordinate of r2 value on figure
        scatter(x,y, 'filled')
        R=corrcoef(x,y);
        R_squared=R(2)^2; %Calc r2
        %text(x_text, y_text, ['R^2 = ' num2str(R_squared)])
        if j==1,
            %title({model_names{i}, p_names{j}})
            title({model_names{i}, p_names{j}, ['R^2 = ' num2str(R_squared)]})
            xlabel('Orginal Parameters')
            ylabel('Fitted Parameters')
        else
            %title(p_names{j});
            title({p_names{j}, ['R^2 = ' num2str(R_squared)]});
        end
        ax2 = subplot(n_params,2,j*2);
        %plot(1:n_subjs,(org_params(:,j)-fitted_params(:,j)),'-o','MarkerSize',3.5,'MarkerFaceColor','b')
        residual=x-y;
        c = residual;
        c(-std(residual)<residual<std(residual))=0;
        scatter(x,residual,[],residual, 'filled')
%         colormap(ax2,'winter');
        colorbar;
        if j==1,
            title({'Residuals', p_names{j}})
            xlabel('Orginal Parameters')
            ylabel('Orginal - Fixed')
        else
            title(p_names{j});
        end
        %if j==1, title('Residuals'); end
    end
    
    fig_name = [model_names{i} '_UniformScatterResidual'];
    %Save the figure
    save_fig(h,fig_name,'fig');
    
end

