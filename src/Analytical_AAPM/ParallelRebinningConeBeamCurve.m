% Rebin the spiral cone-beam projections into cone-parallel geometry
function Projection = ParallelRebinningConeBeamCurve(Projection, useGPU)

NumberOfSet = Projection.NumberOfSet;
for setIdx = 1 : NumberOfSet
    SO = Projection.DataSet{setIdx}.rho(1,1);
    SD = Projection.conf.ConstantRadialDistance;
    YL = Projection.conf.NumberofDetectorColumns;
    ZL = Projection.conf.NumberofDetectorRows;
    dYA = atan(Projection.conf.DetectorElementTransverseSpacing / 2 / SD) * 2;
    YLC = Projection.conf.DetectorCentralElement(1,1);
    
    Proj = Projection.DataSet{setIdx}.projection;
    Proj = permute(Proj,[2 1 3]);
    
    ViewN = size(Projection.DataSet{setIdx}.projection, 3);
    DeltaFai = 2 * pi / Projection.DataSet{setIdx}.NumberofSourceAngularSteps;
    DeltaTheta = DeltaFai;
    DeltaT = Projection.DataSet{setIdx}.deltaT;%tan(Projection.conf.FanAngle * 0.5) * SO * 2 / YL;% 0.9437;%tan(dYA * YL / 2) * SD * 2 / YL;
    PLC = YLC;
    PProj = ParallelRebinCBCurve(single(Proj), int32(YL), int32(ZL),...
       int32(ViewN), single(YLC), single(dYA),...
       single(DeltaTheta), single(PLC), single(DeltaT),...
       single(DeltaFai), single(SO),int32(useGPU), int32(16), int32(32), int32(2));    
    PProj = permute(double(PProj),[2 1 3]);
    Projection.DataSet{setIdx}.projection = PProj;
end
