function [outProj, detWidth] = rebinToPanel(proj, detCenterIdx, ... Center Index (0 based)
    detCellWidth, sourToDet)
[detNum, detHei, prjNum] = size(proj);
% Get Larger part
oneFanAngle = atan(detCellWidth / 2 / sourToDet) * 2.0;
largerAngle = max((detCenterIdx + 0.5) * oneFanAngle, (detNum - 1 - detCenterIdx + 0.5) * oneFanAngle);
%coverAngle = largerAngle * 2;
detWidth = tan(largerAngle) * sourToDet * 2;
newDetCellWidth = detWidth / detNum;
outProj = zeros(size(proj));

for ii = 1 : detNum
    % Fan Angle 
    fanAngle = atan((ii - (detNum + 1) / 2) * newDetCellWidth / sourToDet);
    idx = fanAngle / oneFanAngle + detCenterIdx + 1;
    lowIdx = floor(idx);
    higIdx = ceil(idx);
    if(lowIdx < 1)
        lowIdx = 1;
    end
    if(higIdx < 1)
       continue;
    end
    if(lowIdx > detNum)
        continue;
    end
    if(higIdx > detNum)
        higIdx = detNum;
    end
    lowVal = proj(lowIdx,:,:);
    higVal = proj(higIdx,:,:);
    
    outProj(ii,:,:) = lowVal * (higIdx - idx) + higVal * (idx - lowIdx);
    
end

end