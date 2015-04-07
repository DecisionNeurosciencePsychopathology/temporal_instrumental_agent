function [cost, rts, w, w_update, u_all] = weightfit(params, nbasis, ntrials, func, tvec, edge_correct, plots, rand_sig)
%accept a single response time (rt)
%update the weights w on the basis of temporal function func
%given a radial basis set defined by a vector of means (mu) and a fixed sd (sigma)

    %number of basis functions will have to be specified in a loop outside of weightfit... Too wonky to try to
    %optimize an integer parameter.
    
    %params(1) = proportion margin offset
    %params(2) = sd separation between basis functions (Cohen's d)
    
    rng(111); %predictable random draws for gaussian noise
    if nargin < 6
        edge_correct = 0; %no edge correction
    end
    
    if nargin < 7
        plots = 0;
    end
    
    if nargin < 8
        rand_sig = 0; %sigma for random Gaussian noise
    end
    
    w = zeros(ntrials, nbasis); %initialize zero sampling weights
    
    if strcmpi(func, 'gauss_auc') || strcmpi(func, 'gauss_mult')
        %gauss auc with truncated basis and spread function is the winner.
        %Here, allow for 4 params for optimization: prop_offset, basis_sep, prop_spread, prop_noise
        %prop_spread is the SD of the gaussian update function relative to the width of the interval (i.e., ranges 0..1)
        %prop_noise is the SD of the gaussian noise added to responses to avoid trivial troughs/peaks in basis relative
        %to width of interval (i.e., ranges 0..1)
        prop_offset=params(1);
        basis_sep=params(2);
        prop_spread=params(3);
        if (length(params) == 4), prop_noise=params(4);
        else prop_noise = .04; end
        %override rand_sig input parameter for gauss_auc for now
        rand_sig=prop_noise*range(tvec);
    else        
        if length(params) == 1
            prop_offset = .125; %fixed offset of 12.5%
            basis_sep=1.5;
            temporal_param=params(1);
        elseif length(params) == 2 %fixed separation of 1.5
            if strcmpi(func, 'gauss_point')
                %only two params for gausspoint: offset and sep
                prop_offset=params(1);
                basis_sep=params(2);
            else
                prop_offset = params(1);
                basis_sep=1.5;
                temporal_param=params(2);
            end
        elseif length(params) == 3
            prop_offset = params(1);
            basis_sep = params(2);
            temporal_param=params(3);
        end
    end
    

    %tvec is integer-valued vector of timesteps at which to evaluate rbf cost (e.g., 1:500)
    
    margin_offset = (max(tvec) - min(tvec))*prop_offset; %compute lowest basis center
    
    tmin = min(tvec) - margin_offset; tmax=max(tvec) + margin_offset;
    c=tmin:(tmax-tmin)/(nbasis-1):tmax; %means (centers) of basis functions
    sig = (c(2) - c(1))/basis_sep; %standard deviation of basis functions
    
    %disp(prop_offset);
    c_bounded = c;

    %basis function matrix (extend range to lowest and highest bounds for extreme basis functions)
    t_all=round(tmin-4*sig):round(tmax+4*sig);
    gaussmat_extend = zeros(nbasis,length(t_all));
    for j = 1:nbasis
       gaussmat_extend(j,:) = gaussmf(t_all,[sig c(j)]);
    end
    maxauc=sum(gaussmat_extend,2)*ones(1,length(t_all));
    
    gaussmat_extend=gaussmat_extend./maxauc;
    t_all=tvec;
    
    %gaussmat within observed interval for testing truncated gaussian update
    gaussmat = zeros(nbasis,length(tvec));
    for j = 1:nbasis
        gaussmat(j,:) = gaussmf(tvec,[sig c(j)]);
    end
    
    %normalize gauss functions to have AUC = 1.0
    maxauc=max(sum(gaussmat,2));
    gaussmat_auc1=gaussmat./maxauc;
        
    %normalize gauss functions to have AUC = 1.0 within observed time interval
    %this is essentially a truncated Gaussian basis such that AUC = 1.0 for basis functions within interval
    maxauc=sum(gaussmat,2)*ones(1,length(tvec)); %outer product of vectors to allow for col-wise division below
    gaussmat_trunc=gaussmat./maxauc;

    %equalize the AUC across gaussians
    %[~,indmid]=min(abs(c - median(tvec)));
    %trunc_adjust=max(gaussmat(indmid,:))/max(gaussmat_trunc(indmid,:));
    %gaussmat_trunc=gaussmat_trunc*trunc_adjust;
    
    %figure(10); plot(tvec, gaussmat');
    %figure(11); plot(tvec, gaussmat_trunc');
    %figure(12); plot(tvec, gaussmat_auc1');
    %figure(12); plot(t_all, gaussmat_extend');

   
    %compute weight correction vector for radial basis functions to approximate a f(x) === 1.0 flat function 
    %over the finite interval of interest
    init_params=ones(1,length(c));
    upper_bounds=ones(1,length(c)).*50;
    lower_bounds=ones(1,length(c)).*0;
    fmincon_options = optimoptions(@fmincon,'TolX',1e-7,'TolFun',1e-7,'MaxIter',100000, 'MaxFunEvals', 100000); %request a bit higher precision (doesn't really matter here!)
    
    %this version computes weights for regular gaussian basis
    %edge_correct = 3 is the weight correction factor needed to approximate f(x)===1.0 over the observed interval
    if edge_correct==3
        %this version computes weights for truncated gaussian basis
        [weight_correct, cost_flatrbf] = fmincon(@(params) correctrbf(params, 1:500, c, ones(1,length(c)).*sig, [min(tvec) max(tvec)]), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
        %[weight_correct, cost_flatrbf] = fmincon(@(weight_par) correctrbf(weight_par, tvec, c, ones(1,length(c)).*sig, -Inf), init_params, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);
    else
         weight_correct=ones(1,length(c));
    end
    %disp(weight_correct);
    
    
    %interpolate eligibility trace to give weight to basis functions outside of interval 
    %matching the lowest or highest basis inside the interval of interest
    if edge_correct == 1 %flatten generalization function at the function value for min/max of the interval of interest
        c_bounded(c > max(tvec)) = max(tvec);
        c_bounded(c < min(tvec)) = min(tvec);
    elseif edge_correct == 2 %as above, but use update value for lowest and highest basis within the interval
        c_bounded(c > max(tvec)) = c(max(find(c < max(tvec))));
        c_bounded(c < min(tvec)) = c(min(find(c > min(tvec))));
    end
    
    w_update = zeros(ntrials, nbasis);
    rts = zeros(ntrials, 1); %vector of choices
    %u_all = zeros(ntrials, length(tvec));
    u_all = zeros(ntrials, length(t_all));
    %choose in the middle of the interval to start
    rts(1) = (max(tvec) - min(tvec) + 1)/2;
    
    %fprintf('nbasis: %d, rbfsigma: %.3f\n', nbasis, sig);
    %fprintf('basis centers: %s\n', num2str(c_bounded));
    fprintf('params: %s\n', num2str(params));
    
    if strcmpi(func, 'gauss_auc') || strcmpi(func, 'gauss_mult')
        %compute sum under the curve of a non-truncated gaussian for use in normalizing temporal spread
        sig_update = prop_spread*range(tvec); %sigma as a fraction of the range of the interval
        refspread = sum(gaussmf([min(tvec)-range(tvec):max(tvec)+range(tvec)], [sig_update, median(tvec)]));
    end
    
    for i = 1:ntrials
        %compute generalization function
        if strcmpi(func, 'triangle')
            %inefficient to recompute for each loop step...
            decay_zero = temporal_param; %timesteps until function decays to zero
            rt_diffs = abs(rts(i)-c_bounded);
            w_temp = abs(decay_zero - rt_diffs)/decay_zero;
            w_temp(rt_diffs > decay_zero) = 0; %zero out anything beyond the decay
            w_update(i,:) = w_temp;
            %fprintf('decay_zero: %d\n', decay_zero);
        elseif strcmpi(func, 'exponential')
            lambda = temporal_param;
            w_update(i,:) = lambda.^(abs(rts(i)-c_bounded)); %exponential update
        elseif strcmpi(func, 'square')
            %square function within bounds
            decay_zero = temporal_param; %timesteps until function decays to zero
            rt_diffs = abs(rts(i) - c_bounded);
            w_update(i,:) = (decay_zero - rt_diffs) > 0;
        elseif strcmpi(func, 'square_auc')
            %square function where parameter sets boundary for calculating AUC for each basis function
            decay_zero = temporal_param;
            [~,update_min] = min(abs(rts(i)-decay_zero - t_all)); %min absolute deviation from lower bound of square
            [~,update_max] = min(abs(rts(i)+decay_zero - t_all)); %min absolute deviation from lower bound of square
            %compute sum (AUC) for each basis function between update_min and update_max
            w_temp=sum(gaussmat(:,update_min:update_max),2);
            if ismember(edge_correct, [1 2])
                %use weights of minimum and maximum basis functions within observed interval
                w_temp(c < min(tvec),:) = w_temp(min(find(c >= min(tvec))), :);
                w_temp(c > max(tvec),:) = w_temp(max(find(c <= max(tvec))), :);
            end
            w_update(i,:) = w_temp;
        elseif strcmpi(func, 'exp_auc')
            lambda = temporal_param;
            elig = lambda.^(abs(rts(i)-t_all)); %exponential update
            w_temp = sum(repmat(elig,nbasis,1).*gaussmat, 2);

            if ismember(edge_correct, [1 2])
                %use weights of minimum and maximum basis functions within observed interval
                w_temp(c < min(tvec),:) = w_temp(min(find(c >= min(tvec))), :);
                w_temp(c > max(tvec),:) = w_temp(max(find(c <= max(tvec))), :);
            end
            
            w_update(i,:) = w_temp;
        elseif strcmpi(func, 'gauss_mult')
            
            %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
            elig = gaussmf(tvec, [sig_update, rts(i)]);
            %compute sum of area under the curve of the gaussian function
            auc=sum(elig); %*ones(1,length(elig)); %outer product of vectors to allow for col-wise division below
            
            %divide gaussian update function by its sum so that AUC=1.0
            %note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
            %will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
            elig=elig/auc*refspread;
            %elig=elig/auc;
            
            %N.B.: The multiplication step here works better
            
            %to find the area under the curve of the intersection of two Gaussians, need to take the min of each 
            %PDF (AUC = 1.0) and integrate that over the interval. Here, for the moment, both gaussmat_trunc and
            %elig are PDFs (AUC=1.0) and we sum the mins over all timesteps for the comparison of each basis PDF with
            %the elig PDF. This then represents the proportion overlap between the RBF Gaussian and the spread Gaussian,
            %which will necessarily vary from 0..1.
            %w_temp = arrayfun(@(x) sum(min(gaussmat_trunc(x,:), elig)), 1:size(gaussmat_trunc,1));
            
            %w_basis = repmat(elig,nbasis,1).*gaussmat_auc1;
            w_basis = repmat(elig,nbasis,1).*gaussmat_trunc;
            w_temp = sum(w_basis, 2); %sum area under each basis function's eligibility trace
                        
            if ismember(edge_correct, [1 2])
                %use weights of minimum and maximum basis functions within observed interval
                w_temp(c < min(tvec),:) = w_temp(min(find(c >= min(tvec))), :);
                w_temp(c > max(tvec),:) = w_temp(max(find(c <= max(tvec))), :);
            end
            w_update(i,:) = w_temp;
        elseif strcmpi(func, 'gauss_auc')
            %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
            %elig = gaussmf(tvec, [sig_update, rts(i)]);
            elig = gaussmf(t_all, [sig_update, rts(i)]);
            %compute sum of area under the curve of the gaussian function
            auc=sum(elig); %*ones(1,length(elig)); %outer product of vectors to allow for col-wise division below
            
            %divide gaussian update function by its sum so that AUC=1.0
            %note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
            %will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
            elig=elig/auc;
            
            %to find the area under the curve of the intersection of two Gaussians, need to take the min of each 
            %PDF (AUC = 1.0) and integrate that over the interval. Here, for the moment, both gaussmat_trunc and
            %elig are PDFs (AUC=1.0) and we sum the mins over all timesteps for the comparison of each basis PDF with
            %the elig PDF. This then represents the proportion overlap between the RBF Gaussian and the spread Gaussian,
            %which will necessarily vary from 0..1.
            %w_temp = arrayfun(@(x) sum(min(gaussmat_trunc(x,:), elig)), 1:size(gaussmat_trunc,1));
            %w_temp = arrayfun(@(x) sum(min(gaussmat_auc1(x,:), elig)), 1:size(gaussmat_auc1,1));
            
            w_temp = arrayfun(@(x) sum(min(gaussmat_extend(x,:), elig)), 1:size(gaussmat_extend,1));
                        
            if ismember(edge_correct, [1 2])
                %use weights of minimum and maximum basis functions within observed interval
                w_temp(c < min(tvec)) = w_temp(min(find(c >= min(tvec))));
                w_temp(c > max(tvec)) = w_temp(max(find(c <= max(tvec))));
            end
            w_update(i,:) = w_temp;

        elseif strcmpi(func, 'gauss_point')
            w_temp = gaussmat_trunc(:,find(t_all == rts(i)));
            if ismember(edge_correct, [1 2])
                %use weights of minimum and maximum basis functions within observed interval
                w_temp(c < min(tvec),:) = w_temp(min(find(c >= min(tvec))), :);
                w_temp(c > max(tvec),:) = w_temp(max(find(c <= max(tvec))), :);
            end
            w_update(i,:) = w_temp;
        elseif strcmpi(func, 'gauss')
            %gaussian update function
            sig_update = temporal_param;
            w_update(i,:) = gaussmf(c_bounded, [sig_update, rts(i)]);
        else
            error(['unknown update function: ', func]);
        end
        
        % eh(i,:) = 1./(abs(rts(i)-c)+1); %inverse update
        %w(i+1,:) = w(i,:) + w_update(i,:);
        w(i+1,:) = w(i,:) + w_update(i,:).*weight_correct;
        
        
        %consider adding correction factor?
        %uh(i+1,:) = uh(i,:) + 1.0.*eh(i,:)./correction_factor;
        
        %check that weight correction factor yields f(x)===1.0 within tvec
        %fm=rbfeval(t_all, weight_correct, c, sig);
        %figure(11); plot(tvec,fm(ismember(t_all,tvec)))
        
        %choose next rt
        %scurve = rbfeval(tvec, w(i+1,:), c, ones(1,nbasis).*sig); %sampling curve
        %scurve = rbfeval(t_all, w(i+1,:), c, ones(1,nbasis).*sig); %sampling curve using non-truncated basis
        %scurve = sum(w(i+1,:)'*ones(1,length(tvec)).*gaussmat_trunc);
        scurve = sum(w(i+1,:)'*ones(1,length(tvec)).*gaussmat_auc1);
        %scurve = rbfeval(t_all, w(i+1,:), c, ones(1,nbasis).*sig, [min(tvec), max(tvec)]); %sampling curve
        scurve_inbounds = scurve(ismember(t_all, tvec));
        rt_uncertainty = find(scurve_inbounds==min(scurve_inbounds), 1) + round(rand_sig*randn(1,1));
        if rt_uncertainty < min(tvec), rt_uncertainty=min(tvec);
        elseif rt_uncertainty > max(tvec), rt_uncertainty=max(tvec); end
            
        rts(i+1) = rt_uncertainty; %return position of first min (know least about it)
        u_all(i,:) = scurve;
        
        if plots
            figure(1);
            plot(t_all, u_all(i,:), 'r'); hold on;
            line([min(tvec),min(tvec)],[0,max(u_all(i,:))]);
            line([max(tvec),max(tvec)],[0,max(u_all(i,:))]);
            if i > 1, plot(t_all, u_all(i-1,:), 'b'); end
            scatter(rts(i), 0.2*max(u_all(i,:)), 200);
            hold off;
            drawnow update;
            figure(2);
            %plot(c, w_update(i,:));
            plot(tvec, elig);
            
            %figure(3);
            %plot(t_all, w_basis);
            
            figure(4);
            plot(c, w_update(i,:));
            
            pause(0.1); %slow down plotting
        end
    end

    %only use cost of rbf squared errors on the last trial (after function should have approached random uniform)
    %yhat = rbfeval(tvec, w(size(w,1), :), c, sig)';
    
    %cost function should be deviation from random uniform sampling... Wouldn't that effectively be sum of
    %squared deviations from the mean RBF value within the valid time bounds?
    %cost = sum((yhat - mean(yhat)).^2); %this uses flatness of the u function
    
    %unif = makedist('Uniform','Lower',min(tvec),'Upper',max(tvec));
    %we really want flatness of the rt sampling
    %[~,~,cost]=kstest(rts, [sort(rts) unifcdf(sort(rts), min(tvec), max(tvec))]);
    %[~,cost]=adtest(rts, 'Distribution', unif);%[rts unifcdf(rts, min(tvec), max(tvec))]);
    
    %note that this will collapse bins if expected counts are less than 5, which can lead to identical
    %chi-square values when the number of trials is too small.
    
    %compute chi-square goodness of fit to expected rt counts in random uniform bins
    %divide time interval into 12 equal-sized bins and compute observed versus expected counts for each bin.
    [~,~,stats]=chi2gof(rts,'cdf', @(rts)unifcdf(rts,min(tvec),max(tvec)), 'nbins', 12);
    
    %nbins=12;
    %[~,~,stats]=chi2gof(rts,'expected', ones(1,nbins)*length(rts)/nbins, 'nbins', nbins);
    cost=stats.chi2stat;
    %cost = -1*cost; %maximize p value
    disp(cost);
    
    %cost=sum((rts - mean(rts)).^2);
end
