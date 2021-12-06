load('testBackprojectionIssue');
ReconConf.xn = 256;
ReconConf.yn = 256;
ReconConf.zn = 256;
ReconConf.orgIdxX = 128.5;
ReconConf.orgIdxY = 128.5;
ReconConf.orgIdxZ = 1;

ctimage = Backprojection3DWeighting(Projection, ReconConf, int32(1));
imshow(squeeze(ctimage(128,:,:)),[])