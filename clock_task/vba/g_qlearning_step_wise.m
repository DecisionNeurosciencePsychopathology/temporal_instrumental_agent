function [gx] = g_qlearning_step_wise(x,phi,u,inG)
% simulating a softmax policy to determine furture actions for qlearning

beta = exp(phi);

ntimesteps = inG.ntimesteps;
%Grab only the quit values
%v=x(ntimesteps+1:end);
v=reshape(x,ntimesteps,2);
% if isempty(v_quit)
%     v_quit=x(2)*ones(1,ntimesteps);
% else
%     stop=1;
% end

%Make it so we take the softmax of [q_wait q_quit] at each time step, and the resulting vector of p_quit will be the output
v_max=max(v(:,1),v(:,2));
v_diff = bsxfun(@minus,v,v_max); 
p_choice = (exp(v_diff./beta)) ./ repmat((sum(exp(v_diff./beta),2)),1,2); %Divide by temperature


% p_choice = zeros(length(v),2);
% for i = 1:length(v)
%     p_choice(i,:) = (exp((v(i,:)-max(v(i,:)))/beta)) / (sum(exp((v(i,:)-max(v(i,:)))/beta))); %Divide by temperature
% end
p_quit = p_choice(:,2);
gx = p_quit';
%gx=ones(1,ntimesteps)*.025;



% gx = sig( beta*v_quit );
%  dgdx = zeros(size(x,1),1);
%  dgdx(1:ntimesteps) = beta.*gx.*(1-gx);
%  dgdx(ntimesteps+1:end) = -beta.*gx.*(1-gx);
%  dgdx = dgdx*zeros(1,ntimesteps);
%  dgdP = 0;
 %dgdP = [beta.*v_quit.*gx.*(1-gx)];
