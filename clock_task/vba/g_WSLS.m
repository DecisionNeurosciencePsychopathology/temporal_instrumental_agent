function  [ gx ] = g_WSLS(xt,phi,u,inG)

% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - precision: SD of the Gaussian around previous RT for win-stay
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT

ntimesteps = inG.ntimesteps;

beta = exp(phi(1));                
precision = 1./(1+exp(-phi(2)));

if inG.stay
    stay = 1;
else
    stay = u(2)>0;
end
RT_old = u(1);
%gx = zeros(1,ntimesteps);
x = 1:ntimesteps;

if stay
    %gx(RT_old) = 1;
%     index = gx;
%     index(RT_old) = 1;
    gx = pdf('Normal',x,RT_old,precision)';
    %gx = (exp((index-max(index))/beta)) / (sum(exp((index-max(index))/beta)));
else
    gx = 1/ntimesteps.*ones(1,ntimesteps)';
end
% figure(99);
% plot(gx);
% hold on;
% x = 0:1:ntimesteps;
% y = pdf('Normal',x,RT_old,precision);
% plot(y)
% hold off;
