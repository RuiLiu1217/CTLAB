% analytical Reconstruction demo
% author: Rui Liu

mexcuda Backprojection3DWeight.cu
mexcuda -lcufft ReweightingAndFiltering_GPU.cu
mexcuda ParallelRebinCBCurve.cu

readFromFolder = false;
if readFromFolder
    if ispc
        ProjectionDataPath = 'E:\CT reconstruction\AAPM\L067\L067\full_DICOM-CT-PD';
        dictionaryFile =  'E:\CT reconstruction\AAPM\DICOM-CT-PD-dict_v8.txt';
    elseif isunix
        ProjectionDataPath = '/home/liurui/Desktop/AAPM/AAPM_Data/L067/L067/full_DICOM-CT-PD';
        dictionaryFile =  '/home/liurui/Desktop/AAPM/AAPM_Data/GeneralInfo/GeneralInfo/DICOM-CT-PD-dict_v8.txt';
    else
        error("Platform Mac is not supported");
    end
    Projection = readProjectionData(ProjectionDataPath, dictionaryFile);
else
    load('testProjection.mat');
end

xn = 512;
yn = 512;
zn = 512;
orgIdxX = (xn + 1.0) / 2.0; % MATLAB start with 1 therefore the center index should be 256.5
orgIdxY = (yn + 1.0) / 2.0;
orgIdxZ = 1;
ReconConf = genReconConf(Projection, xn, yn, zn, orgIdxX, orgIdxY, orgIdxZ);
useGPU = 1;

[ctimage, time] = analyticalRecon(Projection, ReconConf, useGPU);
