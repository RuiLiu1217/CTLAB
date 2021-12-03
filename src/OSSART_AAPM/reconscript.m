%-------------------------------------------------------------------------
% Date: Dec. 2, 2021
% Routine: reconscript
% Authors:
%
%   Rui Liu
% Organization:
% Wake Forest Health Sciences & University of Massachusetts Lowell
%
% Aim:
%   The reconstruction script initially defined
%
% Input/Output:
%--------------------------------------------------------------------------
clear;
% compileCOVfuncs; %% Compile the CPP file

%% Read the projection  (NOTE: PLEASE CHANGE THE PATH MANUALLY)
% If there is no 'projectionTest.mat' : [proj, cfg] = readProj('E:\CT reconstruction\AAPM\L067\L067\full_DICOM-CT-PD', 'E:\CT reconstruction\AAPM\DICOM-CT-PD-dict_v8.txt');
[proj, cfg] = readProj('E:\CT reconstruction\AAPM\L067\L067\full_DICOM-CT-PD', 'E:\CT reconstruction\AAPM\DICOM-CT-PD-dict_v8.txt');
% load('projectionTest.mat');
%% Read the image information (NOTE: PLEASE CHANGE THE PATH MANUALLY)
cfgRecon = CollectImageCfg('Info.IMA');

%% Define the number of slices
SLN = 512;

%% The reconstruction configuration
conf = ConvertToReconConf(cfg, cfgRecon, SLN);
conf.recon.dx = 1; % It is suggested being 0.674, I just do the test for the feasibility of OS-SART here.
%% Define the Z positions
h = cfg.SpiralPitchFactor * cfg.DetectorElementAxialSpacing * single(cfg.NumberofDetectorRows) * single(cfg.DetectorFocalCenterRadialDistance) / single(cfg.ConstantRadialDistance);
deltaZ = h / single(cfg.NumberofSourceAngularSteps);

TotalView = cfg.NumOfDataViews;
zPos = linspace(h, deltaZ * TotalView - h, SLN);

%% Rebinning the projection with baodong's parameter
[Proj, Views] = HelicalToFan_routine(proj, cfg, zPos);

%% mask
mask = zeros(conf.recon.XN,conf.recon.YN);
for ii = 1 : conf.recon.XN
    for jj = 1 : conf.recon.YN
        if sqrt(((double(ii) - 0.5 - double(conf.recon.XN) / 2) / (double(conf.recon.XN) / 2))^2 +...
                ((double(jj) - 0.5 - double(conf.recon.YN) / 2) / (double(conf.recon.YN) / 2))^2) < 1.3
            mask(ii,jj) = 1;
        end
           
    end
end
mask = uint8(mask);

useFISTA = 1;
initImg = single(zeros(conf.SLN,conf.recon.XN,conf.recon.YN));
% % % numOfOS_series = [40, 20, 10, 8, 6, 4, 1];
% % % numOfIter_series = [1, 4, 4, 4, 5, 6, 7];
% % % tic;
% % % for ii = 1 : 1
% % %     reconImg = OSSART_AAPM(Proj, Views, initImg, conf, numOfIter_series(ii), numOfOS_series(ii), mask, useFISTA);    
% % %     initImg = reconImg;
% % % end
% % % toc;


numOfIter = 16;
numOfOS = 5;
useFISTA = 1;
initImg = single(zeros(conf.SLN,conf.recon.XN,conf.recon.YN));
reconImg = OSSART_AAPM(Proj, Views, initImg, conf, numOfIter, numOfOS, mask, useFISTA);

