function [gx] = g_qlearning(x,phi,u,inG)
% simulating a softmax policy to determine furture actions for qlearning

beta = exp(phi);

ntimesteps = inG.ntimesteps;
%Grab only the quit values
v_quit=x(ntimesteps+1:end);
% if isempty(v_quit)
%     v_quit=x(2)*ones(1,ntimesteps);
% else
%     stop=1;
% end
p_choice = (exp((v_quit-max(v_quit))/beta)) / (sum(exp((v_quit-max(v_quit))/beta))); %Divide by temperature
gx = p_choice';
%gx=ones(1,ntimesteps)*.025;

%Just checking to see of this breaks
if sum(isnan(gx))>0 || sum(isinf(gx))>0
    stop=1;
end




%Old delete later 
% % epsilon = 1./(1+exp(-P(1)));
% % 
% % % initialize Q with values from x, and actions as a vector of zeros
% % ntimesteps = inF.ntimesteps;
% % Q = [x(1:ntimesteps)  x(ntimesteps+1:end)];
% % actions = zeros(ntimesteps,1);
% % for i = 1:length(Q)
% %     [~, actions(i)] = max(Q(i,:)); %qmax is the maximum value of Q, a is its index/position
% %     
% %     randomNormal = rand(1000,1);
% %     if(randomNormal(randi(length(randomNormal),1)) <= epsilon)
% %         if rand>(1-(1/(ntimesteps-pos.row)))
% %             actions(i)=2;
% %         else
% %             actions(i)=1;
% %         end
% %     end
% % end
