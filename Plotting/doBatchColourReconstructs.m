clear all
close all

Root = '/home/omeacock/Documents/SPRruns/Fig1Runs/';
inBranches = {'Patchy','Nonpatchy','FrozenPatchy','FrozenNonpatchy'};
inTwig = 'SimulationResults.mat';
outBranch = 'ColourCells';
outTwig = 'Frame_%04d.tif';

dispSettings.ImgPath = outTwig;
dispSettings.colourCells = 'Hits';

for r = 1:size(inBranches,2)
    load(fullfile(Root,inBranches{r},inTwig),'endField','trackableData')
    dispSettings.imagedirectory = [Root,inBranches{r},filesep,outBranch];

    paintFieldReconstruction(trackableData,dispSettings,endField)
end