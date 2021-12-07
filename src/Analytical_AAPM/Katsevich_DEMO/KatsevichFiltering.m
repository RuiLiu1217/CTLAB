function FProj=KatsevichFiltering(Proj,ProjScale,DecWidth,DecHeigh,ScanR,StdDis,HelicP,FilterCoeFF,DeltaAngle,FilteringMode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CT/Micro CT lab 
%     Department of Radiology
%     University of Iowa
%  Version of 2003.05.15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Filtering step of katsevich algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ScanR =75;      % (cm) Scanning radius
%StdDis=150;     % (cm) Source to detector distance
%HelicP=25;      % (cm) Helical pitch
%ObjR  =25;      % (cm) Object radius
%DecWidth =107.8;% (cm) Width of detector array
%DecHeigh =39.1; % (cm) Heigh of detector array
%ProjScale=128;  % Number of projection per turn
%RecMatrix = 256;  % Size of reconstructed matrix( for all the three dimensions)  
RebinFactor=1.3;% 0.6<RebinFactor<1.5, The ratio between the resample number and original number in column when filtering 
ParaCoef =FilterCoeFF;% 0<FilterCoef<1, Defining the linear function for constructing the filtering plane
WindowType =1;  % WindowType=1;  rectangle window 
                % WindowType=2;  kaiser window    
                % WindowType=3;  hamming window   
                % WindowType=4;  hanning window   
                % WindowType=5;  blackman window  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate some system parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ProjRes;              % Load  the orginal projection data. 
PI=3.14159265358979;
[YL, ZL, ProjNumber] = size(Proj); % The first dimension is projection row, 
                                   % the second dimension is projecction 
                                   % column and the third is projection views
DeltaL = 2*PI/ProjScale;   
DeltaU = DecWidth/YL;
DeltaV = DecHeigh/ZL;      
HalfZ  = (ZL+1)/2;
HalfY  = (YL+1)/2;

y = ([1:YL]-HalfY)*DeltaU;
z = ([1:ZL]-HalfZ)*DeltaV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Compute the derivative of projections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Compute the derivation of lanbuda%%%
GF = Proj;
for ProjIndex=2:ProjNumber-1
    GF(:,:,ProjIndex) =(Proj(:,:,ProjIndex+1)-...
                        Proj(:,:,ProjIndex-1))/(2*DeltaL);
end
% Deal with the edge problem
  GF(:,:,1) = GF(:,:,2);
  GF(:,:,ProjNumber) = GF(:,:,ProjNumber-1);

 %%%% Compute the derivation of u  %%%%% 
 for Yindex=3:YL-2
    coef = (StdDis^2+(y(Yindex))^2)/(2*DeltaU*StdDis);
    GF(Yindex,:,:) = (GF(Yindex,:,:)+ ...
                      (Proj(Yindex+1,:,:)-Proj(Yindex-1,:,:))*coef);
 end
 % Deal with the edge problem
  temp = Proj(3,:,:)-Proj(1,:,:);
  for Yindex=[1,2]
    coef = (StdDis^2+(y(Yindex))^2)/(2*DeltaU*StdDis);
    GF(Yindex,:,:) = GF(Yindex,:,:)+temp*coef;
  end
  temp = Proj(YL,:,:)-Proj(YL-2,:,:);
  for Yindex=[YL-1,YL]
    coef = (StdDis^2+(y(Yindex))^2)/(2*DeltaU*StdDis);
    GF(Yindex,:,:) = GF(Yindex,:,:)+temp*coef;
  end

%%%%%%% Compute the derivative of v  %%%%%%%%%%%%%
  Y = ndgrid(y,1:ProjNumber);
  Y = reshape(Y,[YL 1 ProjNumber]);
  coef = 1.0/(2*DeltaV*StdDis);
  for Zindex=3:ZL-2
      temp = Y.*(Proj(:,Zindex+1,:)-Proj(:,Zindex-1,:));
      GF(:,Zindex,:) = GF(:,Zindex,:)+z(Zindex).*temp.*coef;
  end
  % Deal with the edge problem
   temp = Y.*(Proj(:,3,:)-Proj(:,1,:));
   GF(:,1,:) = GF(:,1,:)+z(1)*coef*temp;
   GF(:,2,:) = GF(:,2,:)+z(2)*coef*temp;
   temp = Y.*(Proj(:,ZL,:)-Proj(:,ZL-2,:));
   GF(:,ZL-1,:) = GF(:,ZL-1,:)+z(ZL-1)*coef.*temp;
   GF(:,ZL,:)   = GF(:,ZL,:)+z(ZL)*coef.*temp; 
  % disp('Compute the derivative');
  %clear Proj;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Step 2: Resample from the projections and weighted before the 
  %         filtered procedure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Computer the resampling parameter matrix
  AngleRange = 2*PI-DeltaAngle;
  QQ=ceil(ZL*RebinFactor*0.5); 
  AL=2*QQ+1;
  Y = y;
  
  Nata = [1:QQ].'*AngleRange/(QQ-1);
  % angle s2  
  Nata = [-Nata(end:-1:1);0;Nata];
  FilterCoef = ones(QQ,1)*ParaCoef;
  if(FilteringMode==2) 
    FilterCoef =[FilterCoef;0;FilterCoef]; 
  else
    FilterCoef =[1-FilterCoef;0;FilterCoef];  
  end
  AA   = sin(FilterCoef.*Nata)-FilterCoef.*sin(Nata);
  BB   = FilterCoef-1+cos(FilterCoef.*Nata)-FilterCoef.*cos(Nata);
  CC   = sin((1-FilterCoef).*Nata)-sin(Nata)+sin(FilterCoef.*Nata);
  CC(QQ+1)=1;
  % inter distance
  InDs = Nata.*HelicP.*StdDis.*AA./(2*PI*ScanR*CC); 
  % slant ratio
  Rato = Nata.*HelicP.*BB./(2*PI*ScanR*CC);
  Rato(QQ+1) = HelicP/(2*PI*ScanR);
  % Resampling coordinate in the height
  AZCor = InDs(:,ones(YL,1))+Rato*Y;

  % Resample from the projections
  % Resampling results
  ConvRes = zeros(ProjNumber,AL,YL); 
  for Aindex=1:AL
    %sprintf('Sampling the Aindex=%d',Aindex)
    for tempindex=1:YL
      tempV=AZCor(Aindex,tempindex)/DeltaV+HalfZ;
      if tempV<1.01
        ConvRes(:,Aindex,tempindex)=GF(tempindex,1,:);
      elseif tempV>ZL
        ConvRes(:,Aindex,tempindex)=GF(tempindex,ZL,:);
      else
        tpUpv=ceil(tempV);
        tpLwv=tpUpv-1;
        ConvRes(:,Aindex,tempindex)= ...
            GF(tempindex,tpUpv,:)*(tempV-tpLwv)+ ...
            GF(tempindex,tpLwv,:)*(tpUpv-tempV);        
      end %if(tempV<1)
    end %for tempindex=1:YL
  end %for Aindex
  
  % Weighting before filtering
  coef = StdDis./sqrt(StdDis^2+Y(ones(1,AL),:).^2+AZCor.^2);
  for ProjIndex=1:ProjNumber
    ConvRes(ProjIndex,:,:) = squeeze(ConvRes(ProjIndex,:,:)).*coef;
  end
  %disp('Weighting before filtering');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Step 3: Hilbert transform and weighted after the filtered procedure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  nn  = 2^(nextpow2(YL)+1);
  nn2 =2*nn;
  HS = CreateHSP(nn,WindowType);
  tempHS = HS;
  HS=zeros(nn2,1);
  HS(1:nn/2) = tempHS(nn/2+1:nn);
  HS(nn+nn/2+1:nn2) = tempHS(1:nn/2);
  FFT_F = fft(HS);

  % Hilbert filtering
  ConvRes = reshape(ConvRes,ProjNumber*AL,YL);
  for Aindex=1:size(ConvRes,1)
    TempData = double(ConvRes(Aindex,:));
    FFT_S = fft(TempData(:),nn2);
    TempData = real(ifft(FFT_S.*FFT_F));
    ConvRes(Aindex,:) = -TempData(1:YL).';
  end
  ConvRes = reshape(ConvRes,[ProjNumber,AL,YL]);

  % Weighting after filtering
  for ProjIndex=1:ProjNumber
    ConvRes(ProjIndex,:,:) = squeeze(ConvRes(ProjIndex,:,:))./coef;
  end
  %disp('Weighting after filtering');


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Step 4: Rebinning the filtered results
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  FProj = GF;
  for Yindex=1:YL
    %sprintf('Current convert the line Yindex=%d',Yindex)
    % find the minmum and maxmum value of z and corresponding index
    zmax=-DecHeigh;
    zmin= DecHeigh;
    MAXI= 1;
    MINI= AL;
    for tpindex=1:AL
      if(AZCor(tpindex,Yindex)>zmax)
        zmax=AZCor(tpindex,Yindex);
        MAXI=tpindex;
      end
      if(AZCor(tpindex,Yindex)<zmin)
        zmin=AZCor(tpindex,Yindex);
        MINI=tpindex;
      end
    end %for tpindex=1:AL
        
    % rebinning
    BeginAL=0;
    EndAL  =0;
    BeginCoef=0;
    EndCoef  =0;
    for Zindex=1:ZL
        CurrZ =(Zindex-HalfZ)*DeltaV;
      if (CurrZ<zmin)
        BeginAL=MINI;
        EndAL  =MINI;
        BeginCoef=0.5;
        EndCoef  =0.5;
      elseif (CurrZ>zmax)
        BeginAL=MAXI;
        EndAL  =MAXI;
        BeginCoef=0.5;
        EndCoef  =0.5;
      else
        deltamaxB=DecHeigh;
        deltamaxE=DecHeigh;
        for tpindex=1:AL
          if AZCor(tpindex,Yindex)>CurrZ && ...
                AZCor(tpindex,Yindex)-CurrZ<deltamaxE
            deltamaxE=AZCor(tpindex,Yindex)-CurrZ;
            EndAL=tpindex;
          end
          if AZCor(tpindex,Yindex)<CurrZ && ...
                CurrZ-AZCor(tpindex,Yindex)<deltamaxB
            deltamaxB=CurrZ-AZCor(tpindex,Yindex);
            BeginAL=tpindex;
          end
        end %for tpindex=1:AL
        
        if(~((EndAL==BeginAL+1) || (EndAL==BeginAL)))
            %sprintf('Before Adjust QQ=%d, BeginAL=%d, EndAL=%d\n', QQ,BeginAL,EndAL)
            Sign_Dis=CurrZ-Rato(QQ+1)*y(Yindex);            
            if(Sign_Dis>0)
                if(BeginAL<QQ+1)
                     BeginAL=EndAL-1;
                     deltamaxB=CurrZ-AZCor(BeginAL,Yindex);  
                 elseif (EndAL<QQ+1)
                    EndAL=BeginAL+1;
                    deltamaxE=AZCor(EndAL,Yindex)-CurrZ;
                 elseif (EndAL<BeginAL)
                    BeginAL=EndAL-1;
                    deltamaxB=CurrZ-AZCor(BeginAL,Yindex);                    
                 else
                    EndAL=BeginAL+1;
                    deltamaxE=AZCor(EndAL,Yindex)-CurrZ;
                end
            else
                if(BeginAL>QQ+1)
                     BeginAL=EndAL-1;
                     deltamaxB=CurrZ-AZCor(BeginAL,Yindex);  
                 elseif (EndAL>QQ+1)
                     EndAL=BeginAL+1;
                     deltamaxE=AZCor(EndAL,Yindex)-CurrZ;
                 elseif (EndAL>BeginAL)
                    BeginAL=EndAL-1;
                    deltamaxB=CurrZ-AZCor(BeginAL,Yindex);                    
                 else
                    EndAL=BeginAL+1;
                    deltamaxE=AZCor(EndAL,Yindex)-CurrZ;
                end
            end %if(EndAL>QQ+1)  
            %sprintf('After Adjust BeginAL=%d,EndAL=%d\n', BeginAL,EndAL)
        end% if(EndAL~=BeginAL+1)   
        
        BeginCoef= deltamaxE/(deltamaxE+deltamaxB);
        EndCoef  = deltamaxB/(deltamaxE+deltamaxB);
      end %%if CurrZ<Zmin            

      %% rebinning every projections
      FProj(Yindex,Zindex,:) = ConvRes(:,BeginAL,Yindex)*BeginCoef+ ...
          ConvRes(:,EndAL,Yindex)*EndCoef;
    end %for Zindex
  end %for Yindex
  % disp('Rebining after filtering');