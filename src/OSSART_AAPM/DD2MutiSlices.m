%-------------------------------------------------------------------------
% Wake Forest Health Sciences
% Date: Apr, 7, 2016
% Routine: DD2MutiSlices
% Authors:
%
%   Rui Liu (Wake Forest Health)
% Organization:
% Wake Forest Health Sciences
% University of Massachusetts Lowell
%
% Aim:
%   This is a high level Matalb/Freemat wrapper function for the DD forward
%   and back projectors for AAPM competition. The projection is decomposed
%   into the multi slices form.
%
% Input/Output:
%   flag_fw - reference projector / back projector
%       "Proj" / "Back"
%
%
%   cfg
%       contains the geometry parameters of scanning geometry and dimension
%       of the image volume catrecon format is used. The cfg here is
%       different from dd3 cfg in GE. It is reinterpreted and collected
%       from the dicom files in SIEMENS scanning.
%
%
%   input
%       input image or sinogram data. currently only 3D volume is
%       supported. Image or sinogram must have the z-dimension as the first
%       dimension.
%
%
%   view_ind
%       optional, but highly suggested because the projection data in AAPM
%       is too huge to be stored totally in GPU device memory. This allows
%       to specify a subset of views to project.
%
%
%   mask
%       Optional. a 2D byte mask may be provided so that certain region
%       in the image will be ignored from projection
%   startIdx
%       We use four GPUs to generate the projection/backprojection, it is
%       an array indicates the start slice index for each GPU.
function [output] = DD2MutiSlices(flag_fw, cfg, input, hangs, view_ind, mask, startIdx)

XN = int32(cfg.recon.XN);
YN = int32(cfg.recon.YN);
dx = single(cfg.recon.dx);
SLN = int32(cfg.SLN);
sid = single(cfg.acq.sid);

x0 = single(0);
y0 = single(sid);
DNU = int32(cfg.acq.DNU);
PN = int32(cfg.acq.PN);
sdd = single(cfg.acq.sdd);
detCellWidth = single(cfg.acq.detCellWidth);
stepTheta = atan(detCellWidth / 2 / sdd) * 2.0;

detCntIdx = single(cfg.acq.detCntIdx);

imgXCenter = single(0.0);
imgYCenter = single(0.0);

xds = zeros(DNU, 1);
yds = zeros(DNU, 1);
for ii = 1 : DNU
    curBeta = (single(ii) - detCntIdx) * stepTheta; % May be changed, the start index may be 0 or 1.
    xds(ii,1) = sin(curBeta) * sdd;
    yds(ii,1) = sid - cos(curBeta) * sdd;    
end
xds = single(xds);
yds = single(yds);

hangs = single(hangs);

[slnc, viewNum] = size(hangs);

if ~exist('view_ind','var') || isempty(view_ind)
   view_ind = 1 : viewNum; 
end

hangs = single(hangs(:,view_ind));
PN = int32(length(view_ind));

if( ~exist('mask','var') || isempty(mask) )
    mask = int8(ones(XN,YN,2));
else
    mask = int8([mask,mask']);
end

if (~exist('startIdx','var') || isempty(startIdx)) %% According to the configuration in Linux Platform
    numgpu = gpuDeviceCount;
    startIdx = [];
    for ii = 1 : numgpu
        startIdx = [startIdx, int32([ceil((ii - 1) / numgpu * SLN)])];
    end
end


% We use 3 GPUs instead of 4 to calculate because where may be some
% problems if 4 GPUs are all occupied
if strcmp(flag_fw,'Proj')
    output = MultiSlicesProjGPU(x0, y0, DNU, xds, yds, imgXCenter, imgYCenter, ...
        hangs, PN, XN, YN, SLN, single(input), dx, mask, int32(startIdx));
else
    if strcmp(flag_fw,'Back')
        input = input(:,:,view_ind);
        output = MultiSlicesBackGPU(x0, y0, DNU, xds, yds, detCntIdx, imgXCenter, ...
            imgYCenter, hangs, PN, XN, YN, SLN, single(input), dx,...
            mask, int32(startIdx));
    end
end
end