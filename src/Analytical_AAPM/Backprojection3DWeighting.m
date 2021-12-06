function ctimage = Backprojection3DWeighting(Projection, ReconConf, useGPU)

numberOfSet = Projection.NumberOfSet;
XN = ReconConf.xn;
YN = ReconConf.yn;
ZN = ReconConf.zn;

ctimage = zeros(ZN,XN,YN,numberOfSet);

for setIdx = 1 : numberOfSet
    ScanGeom = generateScanGeom(Projection, ReconConf, setIdx);
    Proj = Projection.DataSet{setIdx}.projection;
    ctImg = zeros(ZN, XN, YN);
    Backprojection3DWeight(ScanGeom,Proj,ctImg,int32(useGPU), int32(64), int32(4), int32(1));
    ctimage(:,:,:,setIdx) = ctImg;
end


end




function ScanGeom = generateScanGeom(Projection, ReconConf, setIdx)
    ScanGeom.ScDet(1) = Projection.DataSet{setIdx}.rho(1,1);
    SD = Projection.conf.ConstantRadialDistance;
    ScanGeom.ScDet(2) = SD - ScanGeom.ScDet(1);
    ScanGeom.ScDet(3) = Projection.conf.NumberofDetectorColumns;
    ScanGeom.ScDet(4) = Projection.conf.NumberofDetectorRows;
    ScanGeom.ScDet(5) = Projection.DataSet{setIdx}.deltaT * Projection.conf.NumberofDetectorColumns;
    ScanGeom.ScDet(6) = Projection.conf.DetectorElementAxialSpacing * Projection.conf.NumberofDetectorRows;
    ScanGeom.ScDet(7) = Projection.conf.DetectorCentralElement(1,1);
    ScanGeom.ScDet(8) = Projection.conf.DetectorCentralElement(2,1);
    ScanGeom.ScDet(9) = Projection.conf.SpiralPitchFactor * Projection.DataSet{setIdx}.rho(1,1) / SD;
    %% 
    ScanGeom.Proj(1) = Projection.DataSet{setIdx}.equalentPhi(1);
    ScanGeom.Proj(2) = Projection.DataSet{setIdx}.equalentPhi(end);
    ScanGeom.Proj(3) = size(Projection.DataSet{setIdx}.equalentPhi, 2);
    ScanGeom.Proj(4) = Projection.DataSet{setIdx}.NumberofSourceAngularSteps;
    
    %% 
    ScanGeom.Obj(1) = Projection.conf.DataCollectionDiameter / 2.0;
    ScanGeom.Obj(2) = ReconConf.xn;
    ScanGeom.Obj(3) = ReconConf.yn;
    ScanGeom.Obj(4) = ReconConf.zn;
    ScanGeom.Obj(5) = ReconConf.orgIdxX;
    ScanGeom.Obj(6) = ReconConf.orgIdxY;
    ScanGeom.Obj(7) = ReconConf.orgIdxZ;
    ScanGeom.Obj(8) = ScanGeom.Obj(1) * 2 / ReconConf.xn;
    ScanGeom.Obj(9) = ScanGeom.Obj(1) * 2 / ReconConf.yn;
    ScanGeom.Obj(10) = ScanGeom.Obj(1) * 2 / ReconConf.zn;
    
    %% 
    ScanGeom.Rec(1) = ReconConf.mdelta;
    ScanGeom.Rec(2) = ReconConf.mHSCoef;
    ScanGeom.Rec(3) = ReconConf.mk1;    
    
end


