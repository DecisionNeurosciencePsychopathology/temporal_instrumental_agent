function [rew] = RewFunction(rt, cond) 

Shift = 700;
rt_extended = 7000;

    CEV_frq = 1-((rt+ Shift)./rt_extended);
    CEV_mag= rt_extended.*37./(rt_extended-(rt+Shift));
    IEV_frq= CEV_frq + CEV_frq.*(0.25.*sin((rt*pi)./5000));
    IEV_mag= 2.*CEV_mag -(10.*log(rt+ Shift));
 
if cond ==12 % CEV

    frq = CEV_frq;
    mag = CEV_mag;
    
elseif cond ==34 % CEVR
    
frq=(rt_extended).*37./(200*(rt_extended-(rt+Shift)));
mag=200*(1-((rt+Shift)./ rt_extended));


 elseif cond ==56 % DEV
     
frq= 2*CEV_frq -  IEV_frq;
mag= 10*log(rt+Shift);
  

elseif cond ==78 % IEV
    frq = IEV_frq;
    mag = IEV_mag;
     
end

 if frq>=rand 
     rew =mag; 
 else rew= 0;
 end;

end