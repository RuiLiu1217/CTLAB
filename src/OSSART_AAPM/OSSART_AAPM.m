%-------------------------------------------------------------------------
% Wake Forest Health Sciences
% Date: Apr, 14, 2016
% Routine: OSSART_AAPM
% Authors:
%
%   Rui Liu
% Organization:
% Wake Forest Health Sciences & University of Massachusetts Lowell
%
% Aim:
%   This is the OS-SART function for AAPM with multi-slice fan-beam
%   geometry.
%
% Input/Output:
%   function [Proj, startView] = HelicalToFan(proj, cfg, zPos)
%
%
%   proj - rebined multi slice fan-beam projection data in order SLN,
%   channel direction, projection views
%
%
%   cfg
%       contains the geometry parameters of scanning geometry and dimension
%       of the image volume catrecon format is used. The cfg here is
%       different from dd3 cfg in GE. It is reinterpreted and collected
%       from the dicom files in SIEMENS scanning.
%
%
%   zPos
%       The Z position of the slice we want to get the projection.
%--------------------------------------------------------------------------
function [reconImg] = OSSART_AAPM(Proj, Views, initImg, conf, numOfIter, numOfOS, mask, useFISTA)
[SLN, DNU, numOfViews] = size(Proj);
if( ~exist('useFISTA','var') || isempty(useFISTA))
    useFISTA = 1;
end

if( ~exist('mask','var') || isempty(mask) )
    mask = int8(ones(XN,YN));
else
    mask = int8(mask);
end

if (~exist('numOfOS','var') || isempty(numOfOS))
    numOfOS = 1;
end

view_idx = cell(numOfOS,1);
rowSum = cell(numOfOS,1);
colSum = cell(numOfOS,1);

% Calculate the col sum and row sum;
conf_slice = conf;
conf_slice.SLN = 1;

%% OS-SART weighting value
disp('Calculating the weigthing matrix');
for ii = 1 : numOfOS
    view_idx{ii,1} = (ii:numOfOS:numOfViews)';
    rowSum{ii,1} = DD2MutiSlices('Proj',conf,single(ones(SLN,conf.recon.XN,conf.recon.YN)),Views,view_idx{ii},mask);
    colSum{ii,1} = DD2MutiSlices('Back',conf,single(ones(SLN,conf.acq.DNU,numOfViews)),Views,view_idx{ii},mask);
    rowSum{ii,1} = repmat(rowSum{ii,1}(1,:,:),[SLN,1,1]);
    colSum{ii,1} = repmat(colSum{ii,1}(1,:,:),[SLN,1,1]);
    disp(['Calculating the ', num2str(ii), 'th weighting matrix..']);
end

%% Reconstructed img
reconImg = single(initImg);

% FISTA Weight
t0 = 1;
disp('Reconstrution start');
for ii = 1 : numOfIter
    if useFISTA == 1
        lasImg = reconImg;
    end
    
    for idx = 1 : numOfOS

        nzr = (rowSum{idx,1}~=0);
        nzc = (colSum{idx,1}~=0);
        
        %Project
        prj = DD2MutiSlices('Proj',conf,reconImg,Views,view_idx{idx},mask);
        rep = Proj(:,:,view_idx{idx});
        %Weighting
        prj(nzr) = (rep(nzr) - prj(nzr)) ./ rowSum{idx,1}(nzr);
        
        %Back
        iiddx = 1 : length(view_idx{idx});
        tp = DD2MutiSlices('Back',conf,single(prj),Views(:,view_idx{idx}),iiddx,mask);
        %Weighting
        tp(nzc) = tp(nzc) ./ colSum{idx,1}(nzc);
        
        %
        reconImg = reconImg + tp;
        
        subplot(2,2,1);
        imagesc(squeeze(reconImg(:,256,:)));
        subplot(2,2,2);
        imagesc(squeeze(reconImg(:,:,256)));
        subplot(2,2,3);
        imagesc(squeeze(reconImg(1,:,:)));
        subplot(2,2,4);
        imagesc(squeeze(reconImg(SLN/2,:,:)));
        pause(0.0001);    
        disp(idx);    
    end
    
    if useFISTA == 1
        t1 = (1 + sqrt(1 + 4 * t0 * t0)) / 2;
        reconImg = reconImg + (t0 - 1) / t1 * (reconImg - lasImg);
        t0 = t1;
    end
    
    disp(ii);   
end