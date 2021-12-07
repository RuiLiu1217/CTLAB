%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find Tao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
RatoRr=2/3;
am=2*pi-2*acos(RatoRr);
%am=am*36/35;
%t*cos(am)-cos(t*am)+1-t=0;
tu=0.5;
td=0;
delta=1;
while (delta>0.0001)
    t=(td+tu)/2;
    f=t*cos(am)-cos(t*am)+1-t;
    delta=abs(f);
    if(f>0)
        tu=t;
    else
        td=t;
    end;
end;
t