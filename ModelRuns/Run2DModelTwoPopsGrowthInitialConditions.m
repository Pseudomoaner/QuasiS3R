%Implementation of the model described in Wensink and Lowen 2012: Emergent states in dense systems of active rods: from swarming to turbulence

clear all
close all

homeDirectory = getenv('USERPROFILE');
addpath(genpath([homeDirectory,'\Dropbox\Twork\Code\Bounded 3D Wensink Model']))
addpath(genpath([homeDirectory,'\Dropbox\Twork\Code\export_fig']))
addpath(genpath([homeDirectory,'\Dropbox\Twork\Code\DiffusionTracker']))

%Define file names etc.
RootSim = 'C:\Users\olijm\Desktop\BoundaryTesting';
reconstructionRootName = 'Channel_1';
outputMatName = 'SimulationResults.mat';

%Parameters of model
%Settings applied to entire field - load from file
load([RootSim,'\fieldConfig.mat']);
switch fieldConfig.BoundaryConds %Choose 'periodic' or 'none' (none also covers 'box' from the configuration file)
    case 'none'
        fieldSettings.boundaryConditions = 'none';
    case 'box'
        fieldSettings.boundaryConditions = 'none';
    case 'periodic'
        fieldSettings.boundaryConditions = 'periodic';
end
fieldSettings.fieldWidth = fieldConfig.Width;
fieldSettings.fieldHeight = fieldConfig.Height;
fieldSettings.maxX = fieldConfig.Width;
fieldSettings.maxY = fieldConfig.Height;

fieldSettings.fieldDepth = 10;
fieldSettings.U0 = 250;
fieldSettings.lam = 1.0;
fieldSettings.f0 = 1;
fieldSettings.burnIndt = 0.001;
fieldSettings.burnInSteps = 100;
fieldSettings.colJigRate = 0.01;
fieldSettings.postDivMovement = 'reverse';
fieldSettings.growthRate = 0.015; %Average increase in aspect ratio over one unit of time
fieldSettings.divThresh = 8; %Aspect ratio at which the cell should divide
fieldSettings.zElasticity = inf;
fieldSettings.areaFrac = 0.35;

%Choose your cell and barrier settings
cellSettingsType = 'doubleCell'; %singleCell, doubleCell, LatticedXYCells
barrierSettingsType = 'loaded'; %none, loaded
barrierSettings.CPs = fieldConfig.CPs;
barrierSettings.fieldImg = fieldConfig.FieldDesign;
barrierSettings.resUp = fieldConfig.resUp;

dispSettings.ImgPath1 = 'Frame_Pops_%04d.tif';
dispSettings.ImgPath2 = 'Frame_Orient_%04d.tif';
dispSettings.ImgPath3 = 'Frame_Grey_%04d.tif';
dispSettings.posVec = [100,100,800,800];
dispSettings.cropRect = [0,0,800,800];
dispSettings.colourCells = 'None';

%Settings for the active rods
cellSettings.a1 = 6;
cellSettings.f1 = 1;
cellSettings.r1 = 0;
cellSettings.a2 = 6;
cellSettings.f2 = 0.5;
cellSettings.r2 = 0;

%Output settings
dispSettings.saveFrames = true;

%Processing settings
procSettings.startTime = 1; %Can be useful to cut out some times, if it takes some time for the model to settle down
procSettings.LocalRad = 10;
procSettings.velocitySmoothingSize = 2;
procSettings.minTrackLength = 5;

%Global simulation settings
startMotileDt = 0.1;
samplingRate = 5.0; %How frequently samples of the simulation should be taken
settlingSimTime = 0; %How long it will take for the simulation to settle into an active configuration
targetSimTime = 5000; %Target motile simulation time

%Set up field size (and other variables) for this simulation
startField = WensinkField(fieldSettings.fieldWidth,fieldSettings.fieldHeight,fieldSettings.fieldDepth,fieldSettings.U0,fieldSettings.lam,fieldSettings.boundaryConditions);
startField = startField.populateField(barrierSettingsType,barrierSettings,cellSettingsType,cellSettings,fieldSettings.areaFrac);
% startField = startField.setColours('Position');

procSettings.pixSize = fieldSettings.fieldWidth/(dispSettings.posVec(3)+2); %In the same units as lam

dispSettings.imagedirectory = [RootSim,filesep,'ColourCells'];
if ~exist(dispSettings.imagedirectory,'dir')
    mkdir(dispSettings.imagedirectory);
end

%Make sure that this value of motiledt doesn't make the simulation explode.
tmpFieldSettings = fieldSettings;
tmpFieldSettings.motiledt = startMotileDt;
tmpFieldSettings.FrameSkip = round(samplingRate/tmpFieldSettings.motiledt);
tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
[tmpPCs,~] = simulateWensinkFieldInitial(startField,tmpFieldSettings,dispSettings);
simOK = checkTimestepAppropriate(tmpPCs,tmpFieldSettings);
while ~simOK %Repeatedly reduce simulation step size if needed
    tmpFieldSettings.motiledt = tmpFieldSettings.motiledt;
    tmpFieldSettings.FrameSkip = round(samplingRate/tmpFieldSettings.motiledt);
    tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
    [tmpPCs,~] = simulateWensinkFieldInitial(startField,tmpFieldSettings,dispSettings);
    simOK = checkTimestepAppropriate(tmpPCs,tmpFieldSettings);
end
fieldSettings.motiledt = tmpFieldSettings.motiledt;
fieldSettings.FrameSkip = round(samplingRate/fieldSettings.motiledt);

%Do initial simulation to allow system to reach an active configuration
fieldSettings.motileSteps = ceil(settlingSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
[~,intermediateField] = simulateWensinkFieldInitial(startField,fieldSettings,dispSettings);

%Do another (fully sampled) simulation a longer period of time
fieldSettings.motileSteps = ceil(targetSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
[PCs,endField] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

%Derived measures (for track processing)
fieldSettings.dt = fieldSettings.motiledt * fieldSettings.FrameSkip;
fieldSettings.maxF = round(fieldSettings.motileSteps/fieldSettings.FrameSkip);

[data,trackableData,toMappings,fromMappings] = processModelPCs(PCs,procSettings,fieldSettings);

fullMatOut = [thisSimRoot,filesep,outputMatName];
save(fullMatOut)