function s=calcCost(s,agents,seeds,i,m,o,ntrials,nbasis,ntimesteps,param_idx,clock_options,str)
%In an effort to save lines of repeating code this function takes in
%virtually every parameter known to man and runs each model, which will
%then save the outputed cost in the structure s.

                [s.(agents{1}).(str)(i,o),~,~,~,s.(agents{1}).eV(i,o)] = clock_logistic_operator_kalman_optimize(s.(agents{1}).opt_params(param_idx,:),seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 1, m);
                [s.(agents{2}).(str)(i,o),~,~,~,s.(agents{2}).eV(i,o)] = clock_logistic_operator_kalman_optimize([s.(agents{2}).opt_params(param_idx,:) 0 0.2],seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 0, m);
                [s.(agents{3}).(str)(i,o),~,~,~,s.(agents{3}).eV(i,o)] = clock_logistic_operator_kalman_optimize([s.(agents{3}).opt_params(param_idx,:) 0.9 0.2],seeds(o,:),m.name, ntrials, nbasis, ntimesteps, 0, 0, m);
                
                %Specify agent
                clock_options.agent = 'qlearning';
                clock_options.cond = m.name;
                [s.(agents{4}).(str)(i,o),~,~,~,~,s.(agents{4}).eV(i,o)] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(o,1:2),s.(agents{4}).opt_params(param_idx,:));
                clock_options.agent = 'sarsa';
                [s.(agents{5}).(str)(i,o),~,~,~,~,s.(agents{5}).eV(i,o)] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(o,1:2),s.(agents{5}).opt_params(param_idx,:));
                
               
end



%just in case 
%s.(agents{3}).(str)(i,o)= blah