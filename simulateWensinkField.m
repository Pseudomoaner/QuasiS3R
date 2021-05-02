function [PCs,field] = simulateWensinkField(startField,fS,dS)

PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;

%Actual simulation
for i = 1:fS.motileSteps
    fprintf('Frame is %i\n',i)
    field = field.stepModel(fS.motiledt,fS.f0,fS.zElasticity,fS.growthRate,fS.divThresh,fS.postDivMovement,fS.colJigRate,dS.colourCells);

<<<<<<< Updated upstream
=======
    %Simulate firing events
    if rem(i,fS.FireSkip) == 0
        field = field.calculateHits(fS.fireRange,fS.FireSkip*fS.motiledt,fS.hitRateType);
        field = field.killCells(fS.killThresh,fS.killType);
    end
    
    %Output visualisation
>>>>>>> Stashed changes
    if rem(i,fS.FrameSkip) == 0
        if dS.saveFrames
            imPath = sprintf(dS.ImgPath,fC);
            fullImPath = [dS.imagedirectory, filesep, imPath];
            switch dS.saveType
                case 'draw'
                    outImg = field.drawField();
                    pause(0.01)
                    
                    imwrite(outImg,fullImPath)                    
                case 'plot'
                    outAx = field.plotField(dS.posVec);
                    export_fig(fullImPath,'-m2')
                    
                    cla(outAx)
            end
            disp([fullImPath,' saved.'])
        end
        
        fC = fC + 1;
        
        PCs = interfaceModelAndDiffusionTracker(field,PCs,round(i/fS.FrameSkip)+1);
    end
end