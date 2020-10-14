%Implementation of the model described in Meacock et al 2020: Bacteria
%solve the problem of crowding by moving slowly
clear all
close all

%% File control and parameter definitions

%Define file names etc.
RootSim = 'Choose\Output\Directory';
reconstructionRootName = 'Channel_1';
outputMatName = 'SimulationResults.mat';

%Parameters of model
fieldSettings.fieldWidth = 200;
fieldSettings.fieldHeight = 200;
fieldSettings.maxX = fieldSettings.fieldWidth;
fieldSettings.maxY = fieldSettings.fieldHeight;
fieldSettings.fieldDepth = 10;
fieldSettings.U0 = 250;
fieldSettings.lam = 1.0;
fieldSettings.f0 = 1;
fieldSettings.burnIndt = 0.001;
fieldSettings.burnInSteps = 100;
fieldSettings.colJigRate = 0.0;
fieldSettings.postDivMovement = 'reverse';
fieldSettings.growthRate = 0.0; %Average increase in aspect ratio over one unit of time
fieldSettings.divThresh = 8; %Aspect ratio at which the cell should divide
fieldSettings.zElasticity = inf;
fieldSettings.areaFrac = 0.25;
fieldSettings.boundaryConditions = 'periodic';

%Choose your cell and barrier settings
cellSettingsType = 'LatticedXYCells'; %singleCell, doubleCell, LatticedXYCells
barrierSettingsType = 'none'; %none, loaded
barrierSettings = struct(); %Null version of the barrier settings variable so there's something to pass into the initialisation function

dispSettings.ImgPath = 'Frame_%04d.tif';
dispSettings.colourCells = 'None';

%Settings for the active rods
cellSettings.a = 5;
cellSettings.f = 1;
cellSettings.r = 0;

%Output settings
dispSettings.saveFrames = true;

%Processing settings
procSettings.startTime = 1; %Can be useful to cut out some times, if it takes some time for the model to settle down
procSettings.LocalRad = 10; 
procSettings.velocitySmoothingSize = 2; %Width of the smoothing window applied to the data
procSettings.minTrackLength = 3; %Minimum length a track needs to be to be included in downstream analyses

%Global simulation settings
startMotileDt = 0.1;
samplingRate = 5.0; %How frequently samples of the simulation should be taken
settlingSimTime = 0; %How long it will take for the simulation to settle into an active configuration
targetSimTime = 10; %Target motile simulation time

%Set up field size (and other variables) for this simulation
startField = WensinkField(fieldSettings.fieldWidth,fieldSettings.fieldHeight,fieldSettings.fieldDepth,fieldSettings.U0,fieldSettings.lam,fieldSettings.boundaryConditions);
startField = startField.populateField(barrierSettingsType,barrierSettings,cellSettingsType,cellSettings,fieldSettings.areaFrac);
startField = startField.setColours('Position');

procSettings.pixSize = 0.5; %In the same units as lam. Value is defined by the settings in WensinkField.drawField() (could bring those parameters out to here in the future)

dispSettings.imagedirectory = [RootSim,filesep,'ColourCells'];
if ~exist(dispSettings.imagedirectory,'dir')
    mkdir(dispSettings.imagedirectory);
end

%% Part 1: Make sure that the selected value of motiledt doesn't make the simulation explode, or reduce until you reach numerical stability
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

%% Part 2: Do initial simulation to allow system to reach an active configuration
fieldSettings.motileSteps = ceil(settlingSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
[~,intermediateField] = simulateWensinkFieldInitial(startField,fieldSettings,dispSettings);

%% Part 3: Do another (fully sampled) simulation for a longer period of time - only data from this simulation period will be stored
fieldSettings.motileSteps = ceil(targetSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
[PCs,endField] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

%% Part 4: Process data and save simulation results
fieldSettings.dt = fieldSettings.motiledt * fieldSettings.FrameSkip;
fieldSettings.maxF = round(fieldSettings.motileSteps/fieldSettings.FrameSkip);

[data,trackableData,toMappings,fromMappings] = processModelPCs(PCs,procSettings,fieldSettings);

fullMatOut = [RootSim,filesep,outputMatName];
save(fullMatOut)