clear all
close all

Root = 'D:\SegregationPaper\MixtureSims\PackFrac_0_5';
genericBranch = 'cellForces_%s_%s_aspectRatios_%s_%s_Repeat_1';
inTwig = 'SimulationResults.mat';
outTwig = 'Reconstructs';

F1 = 1;
A1 = 4;

Forces = [2];
AspRats = [5];

Branches = cell(size(Forces,2)*size(AspRats,2),1);

count = 1;
for f = 1:size(Forces,2)
    for a = 1:size(AspRats,2)
        Branches{count} = sprintf(genericBranch,num2str(Forces(f)),num2str(F1),num2str(AspRats(a)),num2str(A1));
        count = count+1;
    end
end

for r = 1:size(Branches,1)
    load([Root,filesep,Branches{r},filesep,inTwig],'fieldSettings','dispSettings','trackableData')
    dispSettings.imagedirectory = [Root,filesep,Branches{r},filesep,outTwig];

    paintFieldReconstructionTwoPops(trackableData,dispSettings,fieldSettings)
end