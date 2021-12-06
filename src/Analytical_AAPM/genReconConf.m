% Generate the reconstruction volume configuration
% NOTE: We always assume that the index of the volume and the projection

function [ReconConf] = genReconConf(Projection, xn, yn, zn, orgIdxX, orgIdxY, orgIdxZ, mk1, mdelta, mHSCoef)

if(~exist('mHSCoef','var') || isempty(mHSCoef))
    ReconConf.mHSCoef = 0.68;
else
    ReconConf.mHSCoef = mHSCoef;
end

if(~exist('mdelta','var') || isempty(mdelta))
    ReconConf.mdelta = 60;
else
    ReconConf.mdelta = mdelta;
end

if(~exist('mk1','var') || isempty(mk1))
    ReconConf.mk1 = 5;
else
    ReconConf.mk1 = mk1;
end

ReconConf.xn = xn;
ReconConf.dx = Projection.conf.DataCollectionDiameter / xn;

ReconConf.yn = yn;
ReconConf.dy = Projection.conf.DataCollectionDiameter / yn;

ReconConf.zn = zn;
ReconConf.dz = ReconConf.dx; % Here we apply the iso-tropic pixel

ReconConf.orgIdxX = orgIdxX;
ReconConf.orgIdxY = orgIdxY;
ReconConf.orgIdxZ = orgIdxZ;

end


