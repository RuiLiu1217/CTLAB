function [Proj, Views] = HelicalToFan_routine(proj, cfg, zPos)
% % HelicalToFan(
% % 		float* Proj,  					// rebinned projection; in order (Channel, View Index, Slice Index)  TODO: need permute after call by MATLAB
% % 		float* Views, 					// rebinned views (View Index, Slice Index) TODO: need permute after call by MATLAB
% % 		float* proj,  					// raw projection data In order : (Height Index, Channel Index, View Index(Total View))
% % 		float* zPos,  					// sampling position
% % 		const int SLN,                  // slice number
% % 		const float SD, 				// source to detector distance
% % 		const float SO,					// source to iso-center distance
% % 		const float BVAngle,        	// The begin view
% % 		const int DetWidth,         	// number of detector columns
% % 		const int DetHeight,        	// number of detector rows
% % 		const float PerDetW,        	// Detector cell size along channel direction
% % 		const float PerDetH,        	// Detector cell size along bench moving direction
% % 		const int DefTimes,         	// Number of views per rotation
% % 		const float DetCenterW,     	// Detector Center Index
% % 		const float SpiralPitchFactor 	// Pitch defined in SIEMENS
% % 		)

SD = cfg.ConstantRadialDistance; % Done
SO = cfg.DetectorFocalCenterRadialDistance; % Done
BVAngle = cfg.DetectorFocalCenterAngularPosition / 180.0 * pi; % Done
DetWidth = cfg.NumberofDetectorColumns;  % Done
DetHeight = cfg.NumberofDetectorRows;  % Done

PerDetW = cfg.DetectorElementTransverseSpacing; % Done
PerDetH = cfg.DetectorElementAxialSpacing; % Done


DefTimes = cfg.NumberofSourceAngularSteps; % Done


DetCenterW = cfg.DetectorCentralElement.X;

%h = cfg.SpiralPitchFactor * cfg.DetectorElementAxialSpacing * cfg.NumberofDetectorRows * SO / SD;

SLN = length(zPos);

[Proj, Views] = HelicalToFanFunc_mex(single(proj),single(zPos),...
    int32(SLN),single(SD),single(SO),single(BVAngle),...
    int32(DetWidth),int32(DetHeight),single(PerDetW),single(PerDetH),...
    int32(DefTimes),single(DetCenterW),single(cfg.SpiralPitchFactor));

Proj = permute(Proj, [3, 1, 2]);
Views = permute(Views,[2, 1]);