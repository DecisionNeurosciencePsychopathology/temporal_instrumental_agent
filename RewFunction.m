function [rew] = RewFunction(rt, cond) 

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
    
elseif strcmpi(cond, 'QUADUP')
    frq = 1 - .62.*rt./5000;
    mag = 0.00002*(rt-1800).^2+20;
    
end

 if frq>=rand 
     rew =mag; 
 else rew= 0;
 end;

end