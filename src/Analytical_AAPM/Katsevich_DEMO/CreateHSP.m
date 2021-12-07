%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CT/Micro CT lab 
%     Department of Radiology
%     University of Iowa
%     Version of 2003.03.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the unit impulse response of the hilbert transform
% XS represents the length which has a form of 2^n
% Index represents the window function type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HS]=CreateHS(XS,index)
Length=XS;
HS=ones(Length,1);
Center=(Length)/2+1;
PI=3.14159265358979;
HS(1)=0;
for i=2:Center-1
   HS(i)=2*sin(PI*(i-Center)/2)^2/(PI*(i-Center));
end;
HS(Center)=0;
for i=Center+1:Length
   HS(i)=2*sin(PI*(i-Center)/2)^2/(PI*(i-Center));
end;
%%Even
%Center=(Length+1)/2;
%for i=1:Length
%    HS(i)=2*sin(PI*(i-Center)/2)^2/(PI*(i-Center));
%end;


switch (index)
case 1
   %rectangle window
   Window=ones(Length,1);
case 2
   %kaiser window
   Window=Kaiser(Length,2.5);
case 3
   %hamming window
   Window=hamming(Length);
case 4
   %hanning window
   Window=hann(Length);
case 5
   %blackman window
   Window=blackman(Length);
end;
HS=HS.*Window;
   
   
   

   
   
   