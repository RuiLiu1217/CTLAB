% Wake Forest Health Sciences
% author: Rui Liu
% Date: Dec. 2, 2021
% The routine that we gather the reconstruction information from the
% projection data
% Input:
%   FileName    : The projection data provided by AAPM in 'dicom'.
%
% Output:
%   cfg         : The necessary information for projection/backprojection
% -----------------------------------------------------------------------
function [cfg] = CollectReconCfg(FileName, dictionaryFile)%, sliceNum)

%% Read the information of the data
info = dicominfo(FileName, 'dictionary', dictionaryFile);

%% HU Calibration Factor. A calibration factor \mu' for the conversion of
% measured linear attenuation coefficients \mu to CT number (mm^-1); CT
% number =1000 * (\mu - \mu') / \mu';
cfg.HUCalibrationFactor = info.HUCalibrationFactor; 

%% The scan field of view
cfg.DataCollectionDiameter = info.DataCollectionDiameter;

%% The pitch
cfg.SpiralPitchFactor = info.SpiralPitchFactor; %0.6 in L067 data

%% The number of detector rows
cfg.NumberofDetectorRows = info.NumberofDetectorRows;

%% The number of detector columns
cfg.NumberofDetectorColumns = info.NumberofDetectorColumns;
                          
%% Detector Element Transverse Spacing
cfg.DetectorElementTransverseSpacing = info.DetectorElementTransverseSpacing;

%% Detector Element Axial Spacing
cfg.DetectorElementAxialSpacing = info.DetectorElementAxialSpacing;

%% Detector Shape (char type)
cfg.DetectorType = info.DetectorShape;

%% Detector Focal Center Angular Position, \phi_0 the azimuthal angles of
% the detector's focal center (rad)
cfg.DetectorFocalCenterAngularPosition = info.DetectorFocalCenterAngularPosition;

%% Detector Focal Center Axial Position, z_0 the z location of the detector's focal center
% the in-plane distances
% between the detector's focal center and the isocenter (mm)
cfg.DetectorFocalCenterAxialPosition = info.DetectorFocalCenterAxialPosition;

%% Detector Focal Center Radial Distance, \rho_0, the in plane distances
% between the detector's focal center and the isocenter
cfg.DetectorFocalCenterRadialDistance = info.DetectorFocalCenterRadialDistance;

%% Detector Central Element: (Column X, Row Y), the index of the detector
% element aligning with the isocenter and the detector's focal center
cfg.DetectorCentralElement.X = info.DetectorCentralElement(1);
cfg.DetectorCentralElement.Y = info.DetectorCentralElement(2);


%% Constant Radial Distance, d_0, the distance between the detector's focal
% center and the detector element specified in Tag(7031,1033) (mm)
% The index start with 1, therefore, when we call the GPU function, we 
% should minus 1, here, the number is not minused.
cfg.ConstantRadialDistance = info.ConstantRadialDistance; 

%% Source Angular Position Shift
% \delta\phi, the \phi offset from the focal sport to the detector's focal
% center (rad) % Maybe not used
cfg.SourceAngularPositionShift = info.SourceAngularPositionShift;

%% Source Axial Position shift
% \delte z, the z offset fro mthe focal spot to the detector's focal center
% (mm). MAYBE NOT USED BUT NOT SURE
cfg.SourceAxialPositionShift = info.SourceAxialPositionShift;

%% Source Radial Distance Shift
% \delta\rho, the \rho offset from the focal spot to the detector's focal
% center (mm). MAYBE NOT USED BUT NOT SURE.
cfg.SourceRadialDistanceShift = info.SourceRadialDistanceShift;


%% Number of source angular steps.
% The number of projections per complete rotation
cfg.NumberofSourceAngularSteps = info.NumberofSourceAngularSteps;

%% Photon Statistics
% An array describing the spatial distribution of photons along the
% direction of the detector columns, from column1 to column M ( neglecting
% the variation across detector rows). Each element of the array
% corresponds to a detector column.
cfg.PhotonStatistics = info.PhotonStatistics;
                            

%% GE Pitch definition
cfg.SpiralPitchGE =  cfg.SpiralPitchFactor * cfg.NumberofDetectorRows;

%% Slice Number
%cfg.SliceNum = sliceNum;
