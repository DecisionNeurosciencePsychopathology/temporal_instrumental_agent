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



% gx = sig( beta*v_quit );
%  dgdx = zeros(size(x,1),1);
%  dgdx(1:ntimesteps) = beta.*gx.*(1-gx);
%  dgdx(ntimesteps+1:end) = -beta.*gx.*(1-gx);
%  dgdx = dgdx*zeros(1,ntimesteps);
%  dgdP = 0;
 %dgdP = [beta.*v_quit.*gx.*(1-gx)];
