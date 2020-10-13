function [PCs,field] = simulateWensinkFieldInitial(startField,fS,dS)
%Separate simulation to allow simulation to settle into active (rather than previous passive) configuration without cell tilting taking place.
PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;
for i = 1:fS.motileSteps
    fprintf('Frame is %i\n',i)
    
    field = field.stepModel(fS.motiledt,fS.f0,inf,fS.growthRate,fS.divThresh,fS.postDivMovement,fS.colJigRate,dS.colourCells);
    
    if rem(i,fS.FrameSkip) == 0
        outImg = field.drawField();
        pause(0.01)
        if dS.saveFrames
            imPath = sprintf(dS.ImgPath,fC);
            fullImPath = [dS.imagedirectory, filesep, imPath];
            
            imwrite(outImg,fullImPath)
            
            disp([fullImPath,' saved.'])
            
            cla
        end
 
        fC = fC + 1;
        
        PCs = interfaceModelAndDiffusionTracker(field,PCs,round(i/fS.FrameSkip)+1);
    end
end