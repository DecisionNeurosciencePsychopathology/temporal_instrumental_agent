function [rew, ev, frq, mag] = RewFunction(rt, cond, userngseed, maxrt)

global rew_rng_state;

if nargin < 3
    userngseed=1; %whether to use global rew_rng_state for draws
end

if nargin < 4, maxrt=4000; end %pertains to the linprob variants where the prob levels off the last 250ms.
minrt=500; %early RT where probability levels off for linprob variants

Shift = 700;
rt_extended = 7000;

CEV_frq = 1-((rt+ Shift)./rt_extended);
CEV_mag= rt_extended.*37./(rt_extended-(rt+Shift));
IEV_frq= CEV_frq + CEV_frq.*(0.25.*sin((rt*pi)./5000));
IEV_mag= 2.*CEV_mag -(10.*log(rt+ Shift));

if strcmpi(cond, 'CEV')
    frq = CEV_frq;
    mag = CEV_mag; 
elseif strcmpi(cond, 'CEVR')
    frq=(rt_extended).*37./(200*(rt_extended-(rt+Shift)));
    mag=200*(1-((rt+Shift)./ rt_extended));
elseif strcmpi(cond, 'DEV')
    frq= 2*CEV_frq -  IEV_frq;
    mag= 10*log(rt+Shift);
elseif strcmpi(cond, 'IEV')
    frq = IEV_frq;
    mag = IEV_mag;
elseif strcmpi(cond, 'QUADUPOLD')
    %previous version where probability was ~1.0 for earliest bump
    frq = 1 - .62.*rt./5000;
    mag = 0.00002*(rt-1800).^2+20;
elseif strcmpi(cond, 'QUADUP')
    %new version where probability of early bump is 0.9 and late bump is 0.34
    %ev ranges from 13.87 to 76.43
    frq = .9 - .56.*rt./5000;
    mag = 0.00002*(rt-1800).^2+20;
elseif strcmpi(cond, 'QUADDOWN')
    %cousin of quadup with bump at 3279ms.
    %ev ranges ranges from 14.49 to 77.14 (pretty close)
    frq = 2e-8*(rt-1800).^2 + .12 + 8e-5*rt;
    mag = -0.00002*(rt-1800).^2 + 224.80;
elseif strcmpi(cond, 'IEVLINPROB')
    rtuppershelf=250; %ms prior to maxrt at which probability levels off

    %for 0-500ms, use the min probability of 0.2
    %for 3750-4000ms, use the max probability of 0.8
    if rt < minrt
        rt = minrt;
    elseif rt > (maxrt - rtuppershelf)
        rt = maxrt - rtuppershelf;
    end
    frq = (rt - minrt)/((maxrt - rtuppershelf - minrt)/0.6) + 0.2; %5416.667 is 3750 - 500 / 0.6
    mag = 1; %1 or 0 outcome
    
elseif strcmpi(cond, 'DEVLINPROB')
    rtuppershelf=250; %ms prior to maxrt at which probability levels off
        
    %for 0-500ms, use the max probability of 0.8
    %for 3750-4000ms, use the min probability of 0.2
    if rt < minrt
        rt = minrt;
    elseif rt > (maxrt - rtuppershelf)
        rt = (maxrt - rtuppershelf);
    end
    
    frq = (maxrt - rtuppershelf - rt)/((maxrt - rtuppershelf - minrt)/0.6) + 0.2;
    mag = 1; %1 or 0 outcome
else
    error(['Unknown function: ' cond]);
end

if userngseed, rng(rew_rng_state); end %draw from reward rng chain
%else rng('shuffle'); end %draw from random chain
rprob=rand;
%fprintf('rew func with prob: %.2f\n', rprob);
if frq >= rprob
    rew = mag;
else rew = 0;
end;
if userngseed, rew_rng_state=rng; end %save state after random draw above

ev=frq*mag; %return expected value

end