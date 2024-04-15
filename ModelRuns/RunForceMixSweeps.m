%Implementation of the model described in Meacock et al 2020: Bacteria
%solve the problem of crowding by moving slowly
clear all
close all

%% File control and parameter definitions

%Define file names etc.
RootSim = '/home/omeacock/Documents/SPRruns/Segregation/MotileMixes/';
reconstructionRootName = 'Channel_1';
outputMatName = 'SimulationResults_F1_%f_F2_%f.mat';

%Settings applied to entire field
fieldSettings.xWidth = 200; %Width of the simulated domain (in units of fieldSettings.lam)
fieldSettings.yHeight = 200; %Height of the simulated domain (in units of fieldSettings.lam)
fieldSettings.postDivMovement = 'reverse'; %How the daughter cell should move following cell division. Either 'reverse' (the opposite direction to the mother) or 'same' (the same direction as the mother).
fieldSettings.fireRange = 2; %Range of the CDI system
fieldSettings.hitRateType = 'distributed'; %Whether CDI hits are diluted over all cells ('distributed') or the per-neighbour hit rate is kept the same regardless of the number of contacts ('constant')
fieldSettings.growthRate = 0.0; %Average increase in aspect ratio over one unit of time.
fieldSettings.divThresh = 8; %Aspect ratio at which the cell should divide.
fieldSettings.areaFrac = 0.35; %Fraction of the total area that should be occupied by cells.

%Choose your cell and barrier settings
barrierSettings.type = 'none'; %Type of static barriers that should be present in simulation - either none or loaded

%Settings for the active rods - note that the use of the 'LatticedXYCells'
%option means that all rods are assumed to be identical. Initialization can
%be customized by writing additional code in the WensinkField.populateField
%function.
cellSettings.type = 'LatticedXYCellsTwoPops'; %Type of rod initialization conditions that should be applied - either singleCell, doubleCell or LatticedXYCells
cellSettings.popFrac = 0.5;
cellSettings.a1 = 5; %Aspect ratio of rods (relative to fieldSettings.lam)
cellSettings.r1 = 0; %Reversal rate associated with each rod
cellSettings.c1 = [1,0,0]; %RGB values for the colour you want to make the cells of population 1
cellSettings.fire1 = 0; %Firing rate of CDI system for this population
cellSettings.pop1 = 'a'; %Population label to specify which cells can kill each other.
cellSettings.a2 = 5;
cellSettings.r2 = 0;
cellSettings.c2 = [0,0.5,1];
cellSettings.fire2 = 0;
cellSettings.pop2 = 'b';

%Output settings
dispSettings.saveFrames = false; %Whether or not to save visualisations of each sampled timepoint
dispSettings.ImgPath = 'Frame_%04d.tif'; %Generic name for each output frame (will be fed into sprintf, so use appropriate string formatting)
dispSettings.colourCells = 'None'; %How rods should be recoloured at each sampling point. If set to 'None', will retain any previously set colour.
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
procSettings.tgtDensity = 1e-2; %Normalized displacement space distance threshold for the link assignment stage of the defect tracking.

%Global simulation settings (defined separately from e.g. field settings so
%they can easily applied uniformly during parameter sweeps).
burnInDt = 0.0625;
startMotileDt = 0.1; %Size of the motility timestep (to begin with)
firingDt = 1; %Size of the timestep for calculating firing events
samplingRate = 5.0; %How frequently samples of the simulation should be taken
burnInSimTime = 0;
targetSimTime = 1000; %Target motile simulation time
contactFind = false; %Whether or not to return structures containing instantaneous cell-cell contact data

forceList = [0.75,1,1.25,1.5,1.75,2];

for f1Ind = 1:size(forceList,2)
    
    cellSettings.f1 = forceList(f1Ind);
    for f2Ind = f1Ind:size(forceList,2)
        cellSettings.f2 = forceList(f2Ind);
        
        %% Part 0: Initialize field for this simulation (including burn-in)
        startField = WensinkField(fieldSettings);
        startField = startField.populateField(barrierSettings,cellSettings,fieldSettings.areaFrac);
        
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
        [~,intermediateField] = simulateWensinkFieldAngularEquilibration(startField,fieldSettings,dispSettings);
        
        %% Part 3: Do another (fully sampled) simulation for a longer period of time - only data from this simulation period will be stored
        fieldSettings.motileSteps = ceil(targetSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
        [PCs,endField] = simulateWensinkField(intermediateField,fieldSettings,dispSettings,contactFind);
        
        %% Part 4: Process data and save simulation results
        fieldSettings.dt = fieldSettings.motiledt * fieldSettings.FrameSkip;
        fieldSettings.maxF = round(fieldSettings.motileSteps/fieldSettings.FrameSkip);
        
        [data,trackableData,toMappings,fromMappings] = processModelPCs(PCs,procSettings,fieldSettings);
        
        [negDefCents,negDefOris,posDefCents,posDefOris] = measureDefects(trackableData,endField,procSettings);

        %For this tracking function to work, ensure FAST's tracking, Progress Bar and helperFuncs folders are added
        %to the path.
        [procDefTracks] = trackDefectsFAST(posDefCents,negDefCents,posDefOris,negDefOris,fieldSettings,procSettings,samplingRate);

        fullMatOut = [RootSim,filesep,sprintf(outputMatName,cellSettings.f1,cellSettings.f2)];
        save(fullMatOut,'data','trackableData','toMappings','fromMappings','fieldSettings','cellSettings','procSettings','samplingRate','startMotileDt','procDefTracks','negDefCents','posDefCents','negDefOris','posDefOris','endField')
    end
end
