%Implementation of the model described in Wensink and Lowen 2012: Emergent states in dense systems of active rods: from swarming to turbulence

clear all
close all

%% File control and parameter definitions

%Define file names etc.
RootSim = 'C:\Users\olijm\Desktop\BoundaryTesting';
reconstructionRootName = 'Channel_1';
outputMatName = 'SimulationResults.mat';

%Parameters of model
%Settings applied to entire field
fieldSettings.fieldWidth = 300; %Width of the simulated domain (in units of fieldSettings.lam)
fieldSettings.fieldHeight = 300; %Height of the simulated domain (in units of fieldSettings.lam)
fieldSettings.maxX = fieldSettings.fieldWidth;
fieldSettings.maxY = fieldSettings.fieldHeight;
fieldSettings.fieldDepth = 10; %Depth of the simulated domain. Value is not critical, provided it is somewhat greater than the length of the longest rod in the simulation
fieldSettings.U0 = 250; %Potential amplitude. Value is not critical, provided it is sufficient to prevent rod-rod crossing (U0 = 250 does this)
fieldSettings.lam = 1.0; %Screening length, defining the interaction distance between rods
fieldSettings.f0 = 1; %Stokesian friction coefficient
fieldSettings.colJigRate = 0.0; %How quickly colours should 'jiggle' noisily (in HSV space). Values above 0 can be useful for visualising lineages of dividing cells. 
fieldSettings.postDivMovement = 'reverse'; %How the daughter cell should move following cell division. Either 'reverse' (the opposite direction to the mother) or 'same' (the same direction as the mother).
fieldSettings.fireRange = 2; %Range of the CDI system
fieldSettings.killType = 'husk'; %Whether to remove cells from simulation after death ('lyse') or inactivate them, but leave their bodies ('husk')
fieldSettings.killThresh = 4; %Number of hits needed to kill a cell
fieldSettings.growthRate = 0.0; %Average increase in aspect ratio over one unit of time.
fieldSettings.divThresh = 8; %Aspect ratio at which the cell should divide.
fieldSettings.zElasticity = 0.6; %Elasticity of the overlying substrate. Set to inf if you want to maintain cells in the monolayer.
fieldSettings.areaFrac = 0.25; %Fraction of the total area that should be occupied by cells.
fieldSettings.boundaryConditions = 'periodic';

%Choose your cell and barrier settings
barrierSettingsType = 'none'; %Type of static barriers that should be present in simulation - either none or loaded
barrierSettings = struct(); %Need to create a dummy variable to pass into the initialization function, even if you don't have any barriers in your system

%Settings for the active rods - note that the use of the 'LatticedXYCells'
%option means that all rods are assumed to be identical. Initialization can
%be customized by writing additional code in the WensinkField.populateField
%function.
cellSettingsType = 'LatticedXYCellsTwoPops'; %Type of rod initialization conditions that should be applied - either singleCell, doubleCell or LatticedXYCells
cellSettings.popFrac = 0.5;
cellSettings.a1 = 4; %Aspect ratio of rods (relative to fieldSettings.lam)
cellSettings.a2 = 5;
cellSettings.f1 = 1.5; %Pushing force applied by each rod
cellSettings.f2 = 5;
cellSettings.r1 = 0; %Reversal rate associated with each rod
cellSettings.r2 = 0;
cellSettings.c1 = [1,1,0];
cellSettings.c2 = [0,1,1];
cellSettings.fire1 = 0;
cellSettings.fire2 = 0;
cellSettings.pop1 = 's';
cellSettings.pop2 = 's';

%Output settings
dispSettings.saveFrames = false; %Whether or not to save visualisations of each sampled timepoint
dispSettings.colourCells = 'None'; %How rods should be recoloured at each sampling point. If set to 'None', will retain any previously set colour.

%Processing settings
procSettings.startTime = 1; %Can be useful to cut out some times, if it takes some time for the model to settle down
procSettings.velocitySmoothingSize = 2; %Size of the smoothing window used to smooth rod position data
procSettings.minTrackLength = 5; %Minimum length of a track to be kept following track assembly
procSettings.pixSize = 0.2; %In the same units as lam. Value is defined by the settings in WensinkField.drawField() (could bring those parameters out to here in the future)

%Global simulation settings (defined separately from e.g. field settings so
%they can easily applied uniformly during parameter sweeps).
burnInDt = 0.0625;
startMotileDt = 0.1; %Size of the timestep (to begin with)
firingDt = 1; %Size of the timestep for calculating firing events
samplingRate = 5.0; %How frequently samples of the simulation should be taken
settlingSimTime = 300; %How long it will take for the simulation to settle into an active configuration
targetSimTime = 600; %Target motile simulation time

%% Part 0: Initialize field for this simulation
startField = WensinkField(fieldSettings.fieldWidth,fieldSettings.fieldHeight,fieldSettings.fieldDepth,fieldSettings.U0,fieldSettings.lam,fieldSettings.boundaryConditions);
startField = startField.populateField(barrierSettingsType,barrierSettings,cellSettingsType,cellSettings,fieldSettings.areaFrac);

%% Part 1: Make sure that the selected value of motiledt doesn't make the simulation explode, or reduce until you reach numerical stability

%Make sure that this value of motiledt doesn't make the simulation explode.
tmpFieldSettings = fieldSettings;
tmpFieldSettings.motiledt = startMotileDt;
tmpFieldSettings.FrameSkip = round(samplingRate/tmpFieldSettings.motiledt);
tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
[tmpPCs,~] = simulateWensinkFieldInitial(startField,tmpFieldSettings,dispSettings);
simOK = checkTimestepAppropriate(tmpPCs,tmpFieldSettings);
simOK = true;
while ~simOK %Repeatedly reduce simulation step size if needed
    tmpFieldSettings.motiledt = tmpFieldSettings.motiledt;
    tmpFieldSettings.FrameSkip = round(samplingRate/tmpFieldSettings.motiledt);
    tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
    [tmpPCs,~] = simulateWensinkFieldInitial(startField,tmpFieldSettings,dispSettings);
    simOK = checkTimestepAppropriate(tmpPCs,tmpFieldSettings);
end
fieldSettings.motiledt = tmpFieldSettings.motiledt;
fieldSettings.FrameSkip = round(samplingRate/fieldSettings.motiledt);
fieldSettings.FireSkip = round(firingDt/fieldSettings.motiledt);

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

exportVtkDataTxts(trackableData,RootSim,fieldSettings,false,[false,false]);

fullMatOut = [RootSim,filesep,outputMatName];
save(fullMatOut)