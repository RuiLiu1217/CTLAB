% % % % % % % % % % % %-------------------------------------------------------------------------
% % % % % % % % % % % % Wake Forest Health Sciences
% % % % % % % % % % % % Date: Mar, 29, 2016
% % % % % % % % % % % % Routine: HelicalToFan
% % % % % % % % % % % % Authors:
% % % % % % % % % % % %
% % % % % % % % % % % %   Baodong Liu, Rui Liu
% % % % % % % % % % % % Organization:
% % % % % % % % % % % % Wake Forest Health Sciences & University of Massachusetts Lowell
% % % % % % % % % % % %
% % % % % % % % % % % % Aim:
% % % % % % % % % % % %   The function that rebinning the helical projection to fan beam projection
% % % % % % % % % % % %   and reconstruct it with linear algebra method or with multi-slice
% % % % % % % % % % % %   fan-beam projection/backprojection
% % % % % % % % % % % %
% % % % % % % % % % % % Input/Output:
% % % % % % % % % % % %   function [Proj, startView] = HelicalToFan(proj, cfg, zPos)
% % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %   proj - helical projection data that we want to rebin
% % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %   cfg
% % % % % % % % % % % %       contains the geometry parameters of scanning geometry and dimension
% % % % % % % % % % % %       of the image volume catrecon format is used. The cfg here is
% % % % % % % % % % % %       different from dd3 cfg in GE. It is reinterpreted and collected
% % % % % % % % % % % %       from the dicom files in SIEMENS scanning.
% % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %   zPos
% % % % % % % % % % % %       The Z position of the slice we want to get the projection.
% % % % % % % % % % % %   
% % % % % % % % % % % %
% % % % % % % % % % % %--------------------------------------------------------------------------
% % % % % % % % % % % function [Proj, Views] = HelicalToFan(proj, cfg, zPos)
% % % % % % % % % % % clc;
% % % % % % % % % % % 
% % % % % % % % % % % SD = cfg.ConstantRadialDistance; % Done
% % % % % % % % % % % SO = cfg.DetectorFocalCenterRadialDistance; % Done
% % % % % % % % % % % BVAngle = cfg.DetectorFocalCenterAngularPosition / 180.0 * pi; % Done
% % % % % % % % % % % DetWidth = cfg.NumberofDetectorColumns;  % Done
% % % % % % % % % % % DetHeight = cfg.NumberofDetectorRows;  % Done
% % % % % % % % % % % 
% % % % % % % % % % % PerDetW = cfg.DetectorElementTransverseSpacing; % Done
% % % % % % % % % % % PerDetH = cfg.DetectorElementAxialSpacing; % Done
% % % % % % % % % % % PerDetA = 2 * atan(0.5 * PerDetW / SD);  
% % % % % % % % % % % 
% % % % % % % % % % % DefTimes = cfg.NumberofSourceAngularSteps; % Done
% % % % % % % % % % % delta = 2 * pi / DefTimes; % The view step size.
% % % % % % % % % % % 
% % % % % % % % % % % DetCenterW = cfg.DetectorCentralElement.X;
% % % % % % % % % % % DetCenterH = cfg.DetectorCentralElement.Y;
% % % % % % % % % % % 
% % % % % % % % % % % DetCenterHreal = (DetHeight-1) / 2 * PerDetH; % What is this? 
% % % % % % % % % % % 
% % % % % % % % % % % h = cfg.SpiralPitchFactor * cfg.DetectorElementAxialSpacing * cfg.NumberofDetectorRows * SO / SD;
% % % % % % % % % % % deltaZ = h / DefTimes;
% % % % % % % % % % % 
% % % % % % % % % % % iii = 1;
% % % % % % % % % % % Proj = cell(length(zPos),1);
% % % % % % % % % % % Views = cell(length(zPos),1);
% % % % % % % % % % % 
% % % % % % % % % % % %% =========================================2D fan-beam rebinning
% % % % % % % % % % % for z = zPos            % set the position of required slice, should be in (h, TotalView * deltaZ - h)
% % % % % % % % % % %     disp(iii);
% % % % % % % % % % %     PData=zeros(DefTimes,DetWidth);
% % % % % % % % % % %     mark=zeros(DefTimes,DetWidth);
% % % % % % % % % % %     start=floor((z-h/2)/deltaZ); %the start projection requaired by rebinning, using the adjacent 2pi segment
% % % % % % % % % % %     ViewAngles = zeros(DefTimes);
% % % % % % % % % % %     for i = 1 : start + DefTimes - 1 %TotalView %read projection, skip the former part
% % % % % % % % % % %         data = proj(:,:,i);
% % % % % % % % % % %         
% % % % % % % % % % %         %% ===================
% % % % % % % % % % %         if i >= start  % only use the required projection. for the other data,only read but don't use it
% % % % % % % % % % %             zz = (i-1) * deltaZ;             % z position
% % % % % % % % % % %             lamda=(i-1) * delta+BVAngle;   % the rotation angle of cone-beam
% % % % % % % % % % %             lamdaFan = mod(lamda,2*pi);    % the rotation angle of fan-beam [0,2pi)
% % % % % % % % % % %             
% % % % % % % % % % %             ibeta = i - start + 1;
% % % % % % % % % % %             ViewAngles(ibeta) = lamdaFan;
% % % % % % % % % % %             
% % % % % % % % % % %           %% ==========================
% % % % % % % % % % %             for j=1:DetWidth
% % % % % % % % % % %                 ang=(j-DetCenterW)*PerDetA;
% % % % % % % % % % %                 SM=SO*cos(ang);  % [Noo,1999, Sigle-slice rebinning method for helical cone-beam CT]
% % % % % % % % % % %                 tanAngV=(z-zz)/SM;
% % % % % % % % % % %                 v=SD*tanAngV;    %curve detector
% % % % % % % % % % %                 dv=(DetCenterHreal+v)/PerDetH+1;
% % % % % % % % % % %                 idv=floor(dv);
% % % % % % % % % % %                 t=dv-idv;
% % % % % % % % % % %                 cosAngV=SD/sqrt(SD^2+v^2);
% % % % % % % % % % %                 %-------------------------linear inerpolation
% % % % % % % % % % %                 temp=0;
% % % % % % % % % % %                 tempmark=0;
% % % % % % % % % % %                 if idv>=1&&idv<=DetHeight&&idv+1>=1&&idv+1<=DetHeight
% % % % % % % % % % %                     temp=data(idv,j)*(1-t)+data(idv+1,j)*(t);
% % % % % % % % % % %                     tempmark=1;
% % % % % % % % % % %                 end
% % % % % % % % % % %                 %-----------------------------
% % % % % % % % % % %                 PData(ibeta,j)=cosAngV*temp;
% % % % % % % % % % %                 mark(ibeta,j)=tempmark;
% % % % % % % % % % %             end
% % % % % % % % % % %         end %% end if
% % % % % % % % % % %     end
% % % % % % % % % % %     
% % % % % % % % % % %     % ----------------------------cut the shape to rectangle
% % % % % % % % % % %     numOfab=0;
% % % % % % % % % % %     for i=1:DefTimes
% % % % % % % % % % %         for j=1:DetWidth
% % % % % % % % % % %             if mark(i,j)==0 %if any ==0, then set the whole row to 0
% % % % % % % % % % %                 numOfab=numOfab+1;
% % % % % % % % % % %                 mark(i,:)=0;
% % % % % % % % % % %                 PData(i,:)=0;
% % % % % % % % % % %                 break;
% % % % % % % % % % %             end
% % % % % % % % % % %         end
% % % % % % % % % % %     end
% % % % % % % % % % %     
% % % % % % % % % % %     ViewAngle=360*(DefTimes-numOfab)/DefTimes; % the rotation angle 
% % % % % % % % % % %     disp(['the total rotation angle (degree)= ' int2str(ViewAngle)]);
% % % % % % % % % % %        
% % % % % % % % % % %     %---------------------------------extract nonzero projections
% % % % % % % % % % %     DefTimesNew=DefTimes-numOfab;
% % % % % % % % % % %     PDataNew=zeros(DefTimesNew,DetWidth);
% % % % % % % % % % %     markNew=zeros(DefTimesNew,1);
% % % % % % % % % % %     ViewAnglesNew=zeros(DefTimesNew,1);
% % % % % % % % % % %     id=0;
% % % % % % % % % % %     for i=1:DefTimes
% % % % % % % % % % %         if mark(i,1)==1
% % % % % % % % % % %             id=id+1;
% % % % % % % % % % %             PDataNew(id,:)=PData(i,:);
% % % % % % % % % % %             markNew(id)=1;
% % % % % % % % % % %             ViewAnglesNew(id)=ViewAngles(i);
% % % % % % % % % % %         end
% % % % % % % % % % %     end
% % % % % % % % % % %     %----------------------------------------------save projection data
% % % % % % % % % % %     Views{iii} = ViewAnglesNew;
% % % % % % % % % % %     Proj{iii} = PDataNew;
% % % % % % % % % % %     iii = iii + 1;
% % % % % % % % % % %     
% % % % % % % % % % % end



%% Modification version according to Dr. Yu's suggestion
function [Proj, Views] = HelicalToFan(proj, cfg, zPos)
clc;

SD = cfg.ConstantRadialDistance; % Done
SO = cfg.DetectorFocalCenterRadialDistance; % Done
BVAngle = cfg.DetectorFocalCenterAngularPosition / 180.0 * pi; % May be changed
DetWidth = cfg.NumberofDetectorColumns;  % Done
DetHeight = cfg.NumberofDetectorRows;  % Done

PerDetW = cfg.DetectorElementTransverseSpacing; % Done
PerDetH = cfg.DetectorElementAxialSpacing; % Done
PerDetA = 2 * atan(0.5 * PerDetW / SD);  

DefTimes = cfg.NumberofSourceAngularSteps; % Done
delta = 2 * pi / DefTimes; % The view step size.

DetCenterW = cfg.DetectorCentralElement.X;
DetCenterH = cfg.DetectorCentralElement.Y;

DetCenterHreal = (DetHeight-1) / 2 * PerDetH; % What is this? 

% % % %  This part is calculated according to Baodong's CPP script
% % h = cfg.SpiralPitchFactor * cfg.DetectorElementAxialSpacing * cfg.NumberofDetectorRows * SO / SD;
% % deltaZ = h / DefTimes;

hi = cfg.SpiralPitchFactor;  %Set the helical pitch
h = hi * PerDetH * DetHeight; %convert to real size
deltaZ = h / DefTimes;

iii = 1;
Proj = cell(length(zPos),1);
Views = cell(length(zPos),1);

%% =========================================2D fan-beam rebinning
for z = zPos            % set the position of required slice, should be in (h, TotalView * deltaZ - h)
    disp(iii);
    PData=zeros(DefTimes,DetWidth);
    mark=zeros(DefTimes,DetWidth);
    start=floor((z-h/2)/deltaZ); %the start projection requaired by rebinning, using the adjacent 2pi segment
    ViewAngles = zeros(DefTimes);
    for i = 1 : start + DefTimes - 1 %TotalView %read projection, skip the former part
        data = proj(:,:,i);
        
        %% ===================
        if i >= start  % only use the required projection. for the other data,only read but don't use it
            zz = (i-1) * deltaZ;             % z position
            lamda=(i-1) * delta+BVAngle;   % the rotation angle of cone-beam
            lamdaFan = mod(lamda,2*pi);    % the rotation angle of fan-beam [0,2pi)
            
            ibeta = i - start + 1;
            ViewAngles(ibeta) = lamdaFan;
            
          %% ==========================
            for j=1:DetWidth
                ang=(j-DetCenterW)*PerDetA;
                SM=SO*cos(ang);  % [Noo,1999, Sigle-slice rebinning method for helical cone-beam CT]
                tanAngV=(z-zz)/SM;
                v=SD*tanAngV;    %curve detector
                dv=(DetCenterHreal+v)/PerDetH+1;
                idv=floor(dv);
                t=dv-idv;
                cosAngV=SD/sqrt(SD^2+v^2);
                %-------------------------linear inerpolation
                temp=0;
                tempmark=0;
                if idv>=1&&idv<=DetHeight&&idv+1>=1&&idv+1<=DetHeight
                    temp=data(idv,j)*(1-t)+data(idv+1,j)*(t);
                    tempmark=1;
                end
                %-----------------------------
                PData(ibeta,j)=cosAngV*temp;
                mark(ibeta,j)=tempmark;
            end
        end %% end if
    end
    
    % ----------------------------cut the shape to rectangle
    numOfab=0;
    for i=1:DefTimes
        for j=1:DetWidth
            if mark(i,j)==0 %if any ==0, then set the whole row to 0
                numOfab=numOfab+1;
                mark(i,:)=0;
                PData(i,:)=0;
                break;
            end
        end
    end
    
    ViewAngle=360*(DefTimes-numOfab)/DefTimes; % the rotation angle 
    disp(['the total rotation angle (degree)= ' int2str(ViewAngle)]);
       
    %---------------------------------extract nonzero projections
    DefTimesNew=DefTimes-numOfab;
    PDataNew=zeros(DefTimesNew,DetWidth);
    markNew=zeros(DefTimesNew,1);
    ViewAnglesNew=zeros(DefTimesNew,1);
    id=0;
    for i=1:DefTimes
        if mark(i,1)==1
            id=id+1;
            PDataNew(id,:)=PData(i,:);
            markNew(id)=1;
            ViewAnglesNew(id)=ViewAngles(i);
        end
    end
    %----------------------------------------------save projection data
    Views{iii} = ViewAnglesNew;
    Proj{iii} = PDataNew;
    iii = iii + 1;
    
end

