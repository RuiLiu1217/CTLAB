% Preweight and filtering the projection
function Projection = PreweightAndFiltering(Projection, useGPU)

numberOfSet = Projection.NumberOfSet;
for setIdx = 1 : numberOfSet
    Proj = Projection.DataSet{setIdx}.projection;
    Proj = permute(Proj,[3 1 2]);
    YL = Projection.conf.NumberofDetectorColumns;
    ZL = Projection.conf.NumberofDetectorRows;
    detCntIdxU = Projection.conf.DetectorCentralElement(1,1);
    detCntIdxV = Projection.conf.DetectorCentralElement(2,1);
    
    
    SO = Projection.DataSet{setIdx}.rho(1,1);
    dYL = Projection.DataSet{setIdx}.deltaT;% tan(Projection.conf.FanAngle * 0.5) * Projection.conf.ConstantRadialDistance * 2 / YL;
    dZL = Projection.conf.DetectorElementAxialSpacing;
    ViewN = size(Proj,1);
    if useGPU == 0
        
        for j = 1 : YL
            t = (j - detCntIdxU) * dYL;
            for k = 1 : ZL
                b = (k - detCntIdxV) * dZL;
                Proj(:,j,k) = Proj(:,j,k)*SO*SO/sqrt(SO^4+(SO*b)^2-(b*t)^2);
            end
        end
        
        % Ramp filtering
        n = -YL:YL;
        hsl=-2./((pi^2)*(4*n.*n-1))/dYL;
        NN = 2^ceil(log(YL*3)/log(2));
        HS = zeros(1,NN);
        HS(1:YL+1)= hsl(YL+1:2*YL+1);
        HS(end-YL+1:end)=hsl(1:YL);
        FtS = fft(HS);
        fpwd = zeros(size(Proj));
        
        parfor i=1:ViewN
            for k=1:ZL
                FtP = fft(Proj(i,:,k),NN);
                tep = ifft(FtS.*FtP);
                fpwd(i,:,k) = real(tep(1:YL));
            end
        end
        Projection.DataSet{setIdx}.projection = permute(fpwd,[2 3 1]);        

    else % use GPU
        Proj = permute(Proj,[2 3 1]);
        fpwd = zeros(size(Proj));
        if(size(Proj,3) < 2048)
            fpwd = double(ReweightingAndFiltering_GPU(single(Proj),single(detCntIdxU), single(detCntIdxV), single(dYL), single(dZL), single(SO)));
        else
            perLen = 2048; % Each time we deal with 2048 views for compatibility. It can be modified with larger numbers.
            segNum = ceil(size(Proj,3) / perLen);
            for ii = 1 : segNum - 1
                rag = 1 + (ii-1) * perLen : ii * perLen;
                fpwd(:,:,rag) = double(ReweightingAndFiltering_GPU(single(Proj(:,:,rag)),single(detCntIdxU), single(detCntIdxV), single(dYL), single(dZL), single(SO)));
            end
            rag = 1 + (segNum - 1) * perLen : (size(Proj,3));
            fpwd(:,:,rag) = double(ReweightingAndFiltering_GPU(single(Proj(:,:,rag)),single(detCntIdxU), single(detCntIdxV), single(dYL), single(dZL), single(SO)));
        end
        
        Projection.DataSet{setIdx}.projection = fpwd;
    end    
end