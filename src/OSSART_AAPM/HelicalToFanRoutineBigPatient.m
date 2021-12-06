function [Proj, Views] = HelicalToFanRoutineBigPatient(proj, standcfg, zPos)

SD = standcfg.sdd;
SO = standcfg.sid;
BVAngle = standcfg.startView;
DetWidth = standcfg.DNU;
DetHeight = standcfg.DNV;

PerDetW = standcfg.colSize;
PerDetH = standcfg.rowSize;

DefTimes = standcfg.ViewPerRot;

DetCenterW = standcfg.colOffSet; % Change the modification

SLN = length(zPos);

pitch = standcfg.PITCH * 0.625 * SD / (standcfg.DNV * standcfg.colSize); % Transfer to pitch


[Proj, Views] = HelicalToFanFunc_mex(single(proj),single(zPos),...
    int32(SLN),single(SD),single(SO),single(BVAngle),...
    int32(DetWidth),int32(DetHeight),single(PerDetW),single(PerDetH),...
    int32(DefTimes),single(DetCenterW),single(pitch));

Proj = permute(Proj, [3, 1, 2]);
Views = permute(Views,[2, 1]);