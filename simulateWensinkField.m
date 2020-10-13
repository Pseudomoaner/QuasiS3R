function [PCs,field] = simulateWensinkField(startField,fS,dS)

PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;

%Actual simulation
for i = 1:fS.motileSteps
    fprintf('Frame is %i\n',i)
    field = field.stepModel(fS.motiledt,fS.f0,fS.zElasticity,fS.growthRate,fS.divThresh,fS.postDivMovement,fS.colJigRate,dS.colourCells);

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
    
    if fS.growthRate > 0
        areaFrac = field.getAreaFraction();
        if areaFrac > fS.areaFrac
            fS.growthRate = 0;
        end
    end
end