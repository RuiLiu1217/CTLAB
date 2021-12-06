% Analytical reconstruction
% Input: 
function [ctimage, time] = analyticalRecon(Projection, ReconConf, useGPU)
if gpuDeviceCount==0
    useGPU = 0;
end
tic;
% Rebinning
disp('Rebinning the cone beam projection');
Projection = ParallelRebinningConeBeamCurve(Projection, useGPU);
% Preweight & filtering the conebeam projection
disp('Preweight the projection');
Projection = PreweightAndFiltering(Projection, useGPU);

% Backprojection
disp('Start Backprojection');
ctimage = Backprojection3DWeighting(Projection, ReconConf, useGPU);
time = toc;


end


