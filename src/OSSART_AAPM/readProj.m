%-------------------------------------------------------------------------
% Wake Forest Health Sciences
% Date: Apr. 5, 2016
% Routine: HelicalToFan
% Author:
%   Rui Liu
%
% Organization:
%   Wake Forest Health Sciences & University of Massachusetts Lowell
%
% Aim:
%   Given the projection path, read the dicoms and rearranged them as
%   Z,XY,ANG order, The projection will be flipped in view dimension. The
%   reading data is suitable for SIEMENS trainning data.
%
% Input/Output:
%   function [proj, cfg] = readProj(FilePath)
% Input:
%   FilePath - The location of dicom format projections
%   
% Output: 
%   proj - helical projection data that we want to rebin
%
%
%   cfg
%       contains the geometry parameters of scanning geometry and dimension
%       of the image volume catrecon format is used. The cfg here is
%       different from dd3 cfg in GE. It is reinterpreted and collected
%       from the dicom files in SIEMENS scanning.
% -------------------------------------------------------------------------
function [proj, cfg] = readProj(FilePath, dictionaryFileName)

fileCluster = [FilePath, '/*.dcm'];
PathInfo = dir(fileCluster);

prjNum = length(PathInfo);

cfg = CollectReconCfg([FilePath,'/',PathInfo(1).name], dictionaryFileName);
cfg.NumOfDataViews = prjNum;

proj = zeros(cfg.NumberofDetectorColumns, cfg.NumberofDetectorRows, prjNum);

parfor prjIdx = 1 : prjNum
    name = [FilePath,'/', PathInfo(prjIdx).name];
    pp = double(dicomread(name)); % Not scaled
    proj(:,:,prjIdx) = pp;
    disp(prjIdx);
end

% Projection is permute to stored in a manner that can be rebined.
proj = permute(proj,[2 1 3]);
proj = proj(:,:,end:-1:1); % Because the differences in GE and SIEMENS definition, we change it to GE standard.

end

