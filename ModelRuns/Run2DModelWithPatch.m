%Implementation of the model described in Meacock et al 2020: Bacteria
%solve the problem of crowding by moving slowly
clear all
close all

%% File control and parameter definitions

%Define file names etc.
RootSim = 'C:\Users\olijm\Desktop\BoundaryTesting';
reconstructionRootName = 'Channel_1';
outputMatName = 'SimulationResults.mat';

%Settings applied to entire field
fieldSettings.fieldWidth = 50; %Width of the simulated domain (in units of fieldSettings.lam)
fieldSettings.fieldHeight = 50; %Height of the simulated domain (in units of fieldSettings.lam)
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
fieldSettings.hitRateType = 'distributed'; %Whether CDI hits are diluted over all cells ('distributed') or the per-neighbour hit rate is kept the same regardless of the number of contacts ('constant')
fieldSettings.growthRate = 0.0; %Average increase in aspect ratio over one unit of time.
fieldSettings.divThresh = 8; %Aspect ratio at which the cell should divide.
fieldSettings.zElasticity = inf; %Elasticity of the overlying substrate. Set to inf if you want to maintain cells in the monolayer.
fieldSettings.areaFrac = 0.35; %Fraction of the total area that should be occupied by cells.
fieldSettings.boundaryConditions = 'periodic';

%Choose your cell and barrier settings
barrierSettingsType = 'none'; %Type of static barriers that should be present in simulation - either none or loaded
barrierSettings = struct(); %Need to create a dummy variable to pass into the initialization function, even if you don't have any barriers in your system

%Settings for the active rods - note that the use of the 'LatticedXYCells'
%option means that all rods are assumed to be identical. Initialization can
%be customized by writing additional code in the WensinkField.populateField
%function.
cellSettingsType = 'LatticedXYCells'; %Type of rod initialization conditions that should be applied - either singleCell, doubleCell or LatticedXYCells
cellSettings.a = 4; %Aspect ratio of rods (relative to fieldSettings.lam)
cellSettings.f = 1; %Pushing force applied by each rod
cellSettings.r = 0; %Reversal rate associated with each rod
cellSettings.c = [1,0,0]; %RGB values for the colour you want to make the cells of population 1
cellSettings.fire = 0; %Firing rate of CDI system for this population
cellSettings.pop = 's'; %Population label to specify which cells can kill each other.

%Settings for the patch initialization, which runs after initial active
%configuration has been reached
patchSettings.fraction = 0.2; %Fraction of rods that should belong to the second population (be part of the patch)
patchSettings.tol = 0.01; %Tolerance for how accurate the patch size should match the specified population fraction. Radius accepted if actual fraction is within range patchSettings.fraction +/- patchSettings.fraction*patchSettings.tol.
patchSettings.force = 1; %Force applied by cells within the patch
patchSettings.reversalRate = 0; %Reversal rate of cells in the patch
patchSettings.colour = [0,0.5,1]; %Colour of cells in patch
patchSettings.fireRate = 0.1; %Firing rate of cells in patch
patchSettings.popLabel = 't'; %Population label of cells in patch

%Output settings
dispSettings.saveFrames = true; %Whether or not to save visualisations of each sampled timepoint
dispSettings.ImgPath = 'Frame_%04d.tif'; %Generic name for each output frame (will be fed into sprintf, so use appropriate string formatting)
dispSettings.colourCells = 'Hits'; %How rods should be recoloured at each sampling point. If set to 'None', will retain any previously set colour.
dispSettings.saveType = 'draw'; %Type of method used to visualise rods - either 'plot' or 'draw'. 'plot' will produce and save a Matlab figure, while 'draw' will draw ellipses directly into an image.
dispSettings.posVec = [100,100,round(900/sqrt(2)),900]; %Determines the location of the plotting figure - only needs to be set if dispSettings.saveType == 'plot'.
dispSettings.imagedirectory = [RootSim,filesep,'ColourCells']; %Defines where the output images will be located
if ~exist(dispSettings.imagedirectory,'dir') %Set up visualisation directory
    mkdir(dispSettings.imagedirectory);
end

%Processing settings
procSettings.startTime = 1; %Can be useful to cut out some times, if it takes some time for the model to settle down
procSettings.velocitySmoothingSize = 2; %Size of the smoothing window used to smooth rod position data
procSettings.minTrackLength = 5; %Minimum length of a track to be kept following track assembly
procSettings.pixSize = 0.2; %In the same units as lam. Value is defined by the settings in WensinkField.drawField() (could bring those parameters out to here in the future)

procSettings.tensorSize = 2; %In the same units as lam. Scale of the smoothing Gaussian filter applied when measuring the orientation field of the output images.
procSettings.incProp = 0.8; %Inclusion proportion for the model training stage of the defect tracking
procSettings.tgtDensity = 1e-3; %Normalized displacement space distance threshold for the link assignment stage of the defect tracking.

%Global simulation settings (defined separately from e.g. field settings so
%they can easily applied uniformly during parameter sweeps).
burnInDt = 0.0625;
startMotileDt = 0.1; %Size of the motility timestep (to begin with)
firingDt = 1; %Size of the timestep for calculating firing events
samplingRate = 5.0; %How frequently samples of the simulation should be taken
burnInSimTime = 0;
settlingSimTime = 200; %How long it will take for the simulation to settle into an active configuration
targetSimTime = 1000; %Target motile simulation time

%% Part 0: Initialize field for this simulation (including burn-in)
startField = WensinkField(fieldSettings.fieldWidth,fieldSettings.fieldHeight,fieldSettings.fieldDepth,fieldSettings.U0,fieldSettings.lam,fieldSettings.boundaryConditions);
startField = startField.populateField(barrierSettingsType,barrierSettings,cellSettingsType,cellSettings,fieldSettings.areaFrac);

tmpFieldSettings = fieldSettings;
tmpFieldSettings.burnIndt = burnInDt;
tmpFieldSettings.burnInSteps = round(burnInSimTime/burnInDt);
tmpFieldSettings.FrameSkip = round(samplingRate/burnInDt);
tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
startField = simulateWensinkFieldBurnIn(startField,tmpFieldSettings,dispSettings);

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

%Setup a patch in the centre of the domain containing the second population
intermediateField = makePatch(intermediateField,patchSettings);

%% Part 3: Do another (fully sampled) simulation for a longer period of time - only data from this simulation period will be stored
fieldSettings.motileSteps = ceil(targetSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
[PCs,endField] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

%% Part 4: Process data and save simulation results
fieldSettings.dt = fieldSettings.motiledt * fieldSettings.FrameSkip;
fieldSettings.maxF = round(fieldSettings.motileSteps/fieldSettings.FrameSkip);

[data,trackableData,toMappings,fromMappings] = processModelPCs(PCs,procSettings,fieldSettings);

fullMatOut = [RootSim,filesep,outputMatName];
save(fullMatOut,'data','trackableData','toMappings','fromMappings','fieldSettings','cellSettings','procSettings','samplingRate','startMotileDt')