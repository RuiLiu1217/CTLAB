% -----------------------------------------------------------------------
% Wake Forest Health Sciences
% Author: Rui Liu
% Date: Apr. 7, 2016
% This is the converter to collect the configuration information from the
% function "CollectImageCfg" and "CollectReconCfg" to a configuration
% struct that store the information can be applied for projection and
% backprojection.
% Input:
%   cfg         : the configuration struct storing the scanning information
%   in AAPM dataset
%
%   cfgRecon    : the configuration struct storing the reconstruction
%   information in AAPM training dataset
%
%   SLN         : slice number
% Output:
%   conf        : the configuration with necessary information for
%   reconstruction
% -----------------------------------------------------------------------
function [conf] = ConvertToReconConf(cfg, cfgRecon, SLN)

conf.recon.XN = cfgRecon.Width;
conf.recon.YN = cfgRecon.Height;

conf.SLN = SLN;
conf.acq.sid = cfg.DetectorFocalCenterRadialDistance;
conf.acq.DNU = cfg.NumberofDetectorColumns;
conf.acq.PN = cfg.NumberofSourceAngularSteps;
conf.acq.sdd = cfg.ConstantRadialDistance;
conf.acq.detCellWidth = cfg.DetectorElementTransverseSpacing;
conf.acq.detCntIdx = cfg.DetectorCentralElement.X;
conf.recon.dx = cfgRecon.PixelSize;
end