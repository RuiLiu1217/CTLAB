% Read the projection data and then divide them into several groups
% according to the flying spot mode
function Projection = readProjectionData(rawDICOMFilePath, dictionaryFileName)
%% Assumption: In the same folder, the scanning mode will not be changed
s = dir([rawDICOMFilePath, '/*.dcm']);
allNames = {s.name};
fileNum = max(size(allNames));
% Get the header
str = [char(rawDICOMFilePath),'/',char(allNames(1,1))]; %
header = dicominfo(str, 'dictionary', dictionaryFileName); % All information are stored in the header, many of them share the same value

%% The shared configuration
Projection.conf.FlyingFocalSpotMode = header.FlyingFocalSpotMode;
Projection.conf.HUCalibrationFactor = double(header.HUCalibrationFactor);
Projection.conf.DataCollectionDiameter = double(header.DataCollectionDiameter);
Projection.conf.SpiralPitchFactor = double(header.SpiralPitchFactor);
Projection.conf.ConstantRadialDistance = double(header.ConstantRadialDistance);
Projection.conf.NumberofDetectorRows = double(header.NumberofDetectorRows);
Projection.conf.NumberofDetectorColumns = double(header.NumberofDetectorColumns);
Projection.conf.DetectorElementTransverseSpacing = double(header.DetectorElementTransverseSpacing);
Projection.conf.DetectorElementAxialSpacing = double(header.DetectorElementAxialSpacing);
Projection.conf.ConstantRadialDistance = double(header.ConstantRadialDistance);
Projection.conf.DetectorCentralElement = double(header.DetectorCentralElement);
Projection.conf.FanAngle = atan(Projection.conf.DetectorElementTransverseSpacing / Projection.conf.ConstantRadialDistance / 2) * 2 * Projection.conf.NumberofDetectorColumns;


PhotonStatistics = zeros(header.NumberofDetectorColumns, fileNum);
projection = zeros(header.NumberofDetectorColumns, header.NumberofDetectorRows, fileNum); % Initially defined projection data
rho = zeros(fileNum, 1);
phi = zeros(fileNum, 1); % This is the left-hand based real projection views
z = zeros(fileNum, 1);

deltaRho = zeros(fileNum, 1);
deltaPhi = zeros(fileNum, 1);
deltaZ = zeros(fileNum,1);


FocalSpotPositions = zeros(3, fileNum);

FocalSpotPositionsX = zeros(fileNum,1);
FocalSpotPositionsY = zeros(fileNum,1);
FocalSpotPositionsZ = zeros(fileNum,1);

parfor fileIndex = 1 : fileNum
    str = [char(rawDICOMFilePath),'/',char(allNames(1,fileIndex))];
    if mod(fileIndex,100) == 0
        disp(['Read file index ', num2str(fileIndex)]);
    end
    hd = dicominfo(str, 'dictionary', dictionaryFileName);
    rho(fileIndex, 1) = double(hd.DetectorFocalCenterRadialDistance);
    deltaRho(fileIndex, 1) = double(hd.SourceRadialDistanceShift);
    phi(fileIndex, 1) = double(hd.DetectorFocalCenterAngularPosition);
    deltaPhi(fileIndex, 1) = double(hd.SourceAngularPositionShift);
    z(fileIndex, 1) = double(hd.DetectorFocalCenterAxialPosition);
    deltaZ(fileIndex, 1) = double(hd.SourceAxialPositionShift);
    
    
    FocalSpotPositionsX(fileIndex,1) = -rho(fileIndex,1) * sin(phi(fileIndex,1));
    FocalSpotPositionsY(fileIndex,1) = rho(fileIndex,1) * cos(phi(fileIndex, 1));
    FocalSpotPositionsZ(fileIndex,1) = z(fileIndex, 1);

    RescaleIntercept = double(hd.RescaleIntercept);
    RescaleSlope = double(hd.RescaleSlope);
    
    PhotonStatistics(:,fileIndex) = double(hd.PhotonStatistics);
    
    
    projection(:,:,fileIndex) = double(dicomread(str)) * RescaleSlope + RescaleIntercept;
end

%% Focal Spot positions
FocalSpotPositions(1,:) = FocalSpotPositionsX';
FocalSpotPositions(2,:) = FocalSpotPositionsY';
FocalSpotPositions(3,:) = FocalSpotPositionsZ';

%projection = projection(:,:,:);

projection = projection(:,end:-1:1,:);


switch Projection.conf.FlyingFocalSpotMode
    case 'FFSNONE'
        Projection.NumberOfSet = 1;
        Projection.DataSet{1}.projection = projection;
        Projection.DataSet{1}.PhotonStatistics = PhotonStatistics;
        Projection.DataSet{1}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps);
        Projection.DataSet{1}.rho = rho;
        Projection.DataSet{1}.phi = phi;
        Projection.DataSet{1}.z = z;
        Projection.DataSet{1}.deltaRho = deltaRho;
        Projection.DataSet{1}.deltaPhi = deltaPhi;
        Projection.DataSet{1}.deltaZ = deltaZ;
        
        Projection.DataSet{1}.FocalSpotPositions = FocalSpotPositions;
        startView = Projection.DataSet{1}.phi(1,1);
        Projection.DataSet{1}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{1}.NumberofSourceAngularSteps;
        Projection.DataSet{1}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{1}.phi)) - 1) * stepView, max(size(Projection.DataSet{1}.phi)));
        
        % Reverse the projection
        
        
        % We need one projection data set
    case 'FFSZ'
        Projection.NumberOfSet = 2;
        Projection.DataSet{1}.projection = projection(:,:,1:2:end);
        Projection.DataSet{1}.PhotonStatistics = PhotonStatistics(:,1:2:end);
        Projection.DataSet{1}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 2;
        Projection.DataSet{1}.rho = rho(1:2:end,1);
        Projection.DataSet{1}.phi = phi(1:2:end,1);
        Projection.DataSet{1}.z = z(1:2:end,1);
        Projection.DataSet{1}.deltaRho = deltaRho(1:2:end,1);
        Projection.DataSet{1}.deltaPhi = deltaPhi(1:2:end,1);
        Projection.DataSet{1}.deltaZ = deltaZ(1:2:end,1);
        
        
        Projection.DataSet{1}.FocalSpotPositions = FocalSpotPositions(:,1:2:end);
        startView = Projection.DataSet{1}.phi(1,1);
        Projection.DataSet{1}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{1}.NumberofSourceAngularSteps;
        Projection.DataSet{1}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{1}.phi)) - 1) * stepView, max(size(Projection.DataSet{1}.phi)));
        
        
        
        Projection.DataSet{2}.projection = projection(:,:,2:2:end);
        Projection.DataSet{2}.PhotonStatistics = PhotonStatistics(:,2:2:end);
        Projection.DataSet{2}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 2;
        Projection.DataSet{2}.rho = rho(2:2:end,1);
        Projection.DataSet{2}.phi = phi(2:2:end,1);
        Projection.DataSet{2}.z = z(2:2:end,1);
        Projection.DataSet{2}.deltaRho = deltaRho(2:2:end,1);
        Projection.DataSet{2}.deltaPhi = deltaPhi(2:2:end,1);
        Projection.DataSet{2}.deltaZ = deltaZ(2:2:end,1);
        
        Projection.DataSet{2}.FocalSpotPositions = FocalSpotPositions(:,2:2:end);
        startView = Projection.DataSet{2}.phi(1,1);
        Projection.DataSet{2}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{2}.NumberofSourceAngularSteps;
        Projection.DataSet{2}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{2}.phi)) - 1) * stepView, max(size(Projection.DataSet{2}.phi)));
        
        
        % We need two projection data sets
    case 'FFSXY'
        Projection.NumberOfSet = 2;
        Projection.NumberOfSet = 2;
        Projection.DataSet{1}.projection = projection(:,:,1:2:end);
        Projection.DataSet{1}.PhotonStatistics = PhotonStatistics(:,1:2:end);
        Projection.DataSet{1}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 2;
        Projection.DataSet{1}.rho = rho(1:2:end,1);
        Projection.DataSet{1}.phi = phi(1:2:end,1);
        Projection.DataSet{1}.z = z(1:2:end,1);
        Projection.DataSet{1}.deltaRho = deltaRho(1:2:end,1);
        Projection.DataSet{1}.deltaPhi = deltaPhi(1:2:end,1);
        Projection.DataSet{1}.deltaZ = deltaZ(1:2:end,1);
        
        Projection.DataSet{1}.FocalSpotPositions = FocalSpotPositions(:,1:2:end);
        startView = Projection.DataSet{1}.phi(1,1);
        Projection.DataSet{1}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{1}.NumberofSourceAngularSteps;
        Projection.DataSet{1}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{1}.phi)) - 1) * stepView, max(size(Projection.DataSet{1}.phi)));
        
        
        Projection.DataSet{2}.projection = projection(:,:,2:2:end);
        Projection.DataSet{2}.PhotonStatistics = PhotonStatistics(:,2:2:end);
        Projection.DataSet{2}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 2;
        Projection.DataSet{2}.rho = rho(2:2:end,1);
        Projection.DataSet{2}.phi = phi(2:2:end,1);
        Projection.DataSet{2}.z = z(2:2:end,1);
        Projection.DataSet{2}.deltaRho = deltaRho(2:2:end,1);
        Projection.DataSet{2}.deltaPhi = deltaPhi(2:2:end,1);
        Projection.DataSet{2}.deltaZ = deltaZ(2:2:end,1);
        Projection.DataSet{2}.FocalSpotPositions = FocalSpotPositions(:,2:2:end);
        startView = Projection.DataSet{2}.phi(1,1);
        Projection.DataSet{2}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{2}.NumberofSourceAngularSteps;
        Projection.DataSet{2}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{2}.phi)) - 1) * stepView, max(size(Projection.DataSet{2}.phi)));
        
        % We need two projection data sets
    case 'FFSXYZ'
        Projection.NumberOfSet = 4;
        
        Projection.DataSet{1}.projection = projection(:,:,1:4:end);
        Projection.DataSet{1}.PhotonStatistics = PhotonStatistics(:,1:4:end);
        Projection.DataSet{1}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 4;
        Projection.DataSet{1}.rho = rho(1:4:end,1);
        Projection.DataSet{1}.phi = phi(1:4:end,1);
        Projection.DataSet{1}.z = z(1:4:end,1);
        Projection.DataSet{1}.deltaRho = deltaRho(1:4:end,1);
        Projection.DataSet{1}.deltaPhi = deltaPhi(1:4:end,1);
        Projection.DataSet{1}.deltaZ = deltaZ(1:4:end,1);
        Projection.DataSet{1}.FocalSpotPositions = FocalSpotPositions(:,1:4:end);
        startView = Projection.DataSet{1}.phi(1,1);
        Projection.DataSet{1}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{1}.NumberofSourceAngularSteps;
        Projection.DataSet{1}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{1}.phi)) - 1) * stepView, max(size(Projection.DataSet{1}.phi)));
        
        
        Projection.DataSet{2}.projection = projection(:,:,2:4:end);
        Projection.DataSet{2}.PhotonStatistics = PhotonStatistics(:,2:4:end);
        Projection.DataSet{2}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 4;
        Projection.DataSet{2}.rho = rho(2:4:end,1);
        Projection.DataSet{2}.phi = phi(2:4:end,1);
        Projection.DataSet{2}.z = z(2:4:end,1);
        Projection.DataSet{2}.deltaRho = deltaRho(2:4:end,1);
        Projection.DataSet{2}.deltaPhi = deltaPhi(2:4:end,1);
        Projection.DataSet{2}.deltaZ = deltaZ(2:4:end,1);
        Projection.DataSet{2}.FocalSpotPositions = FocalSpotPositions(:,2:4:end);
        startView = Projection.DataSet{2}.phi(1,1);
        Projection.DataSet{2}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{2}.NumberofSourceAngularSteps;
        Projection.DataSet{2}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{2}.phi)) - 1) * stepView, max(size(Projection.DataSet{2}.phi)));
        
        
        
        Projection.DataSet{3}.projection = projection(:,:,3:4:end);
        Projection.DataSet{3}.PhotonStatistics = PhotonStatistics(:,3:4:end);
        Projection.DataSet{3}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 4;
        Projection.DataSet{3}.rho = rho(3:4:end,1);
        Projection.DataSet{3}.phi = phi(3:4:end,1);
        Projection.DataSet{3}.z = z(3:4:end,1);
        Projection.DataSet{3}.deltaRho = deltaRho(3:4:end,1);
        Projection.DataSet{3}.deltaPhi = deltaPhi(3:4:end,1);
        Projection.DataSet{3}.deltaZ = deltaZ(3:4:end,1);
        Projection.DataSet{3}.FocalSpotPositions = FocalSpotPositions(:,3:4:end);
        startView = Projection.DataSet{3}.phi(1,1);
        Projection.DataSet{3}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{3}.NumberofSourceAngularSteps;
        Projection.DataSet{3}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{3}.phi)) - 1) * stepView, max(size(Projection.DataSet{3}.phi)));
        
        
        
        Projection.DataSet{4}.projection = projection(:,:,4:4:end);
        Projection.DataSet{4}.PhotonStatistics = PhotonStatistics(:,4:4:end);
        Projection.DataSet{4}.NumberofSourceAngularSteps = double(header.NumberofSourceAngularSteps) / 4;
        Projection.DataSet{4}.rho = rho(4:4:end,1);
        Projection.DataSet{4}.phi = phi(4:4:end,1);
        Projection.DataSet{4}.z = z(4:4:end,1);
        Projection.DataSet{4}.deltaRho = deltaRho(4:4:end,1);
        Projection.DataSet{4}.deltaPhi = deltaPhi(4:4:end,1);
        Projection.DataSet{4}.deltaZ = deltaZ(4:4:end,1);
        
        Projection.DataSet{4}.FocalSpotPositions = FocalSpotPositions(:,4:4:end);
        startView = Projection.DataSet{4}.phi(1,1);
        Projection.DataSet{4}.deltaT = sin(Projection.conf.FanAngle / 2) * Projection.DataSet{1}.rho(1,1) * 2 / Projection.conf.NumberofDetectorColumns;
        stepView = 2 * pi / Projection.DataSet{4}.NumberofSourceAngularSteps;
        Projection.DataSet{4}.equalentPhi = linspace(startView, startView + (max(size(Projection.DataSet{4}.phi)) - 1) * stepView, max(size(Projection.DataSet{4}.phi)));
        
end

end