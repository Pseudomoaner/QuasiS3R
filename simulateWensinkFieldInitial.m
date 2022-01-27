function [PCs,field] = simulateWensinkFieldInitial(startField,fS,dS)
%Separate simulation to allow simulation to settle into active (rather than previous passive) configuration without cell tilting taking place.
PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;
for i = 1:fS.motileSteps
    fprintf('Initialize frame is %i of %i\n',i,fS.motileSteps)
    
    field = field.stepModel(fS.motiledt,fS.growthRate,fS.divThresh,fS.postDivMovement,dS.colourCells,'Initial');
    
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
                    outAx = field.plotField(dS.posVec,dS.showContacts);
                    export_fig(fullImPath,'-m2')
                    
                    cla(outAx)
            end
            disp([fullImPath,' saved.'])
        end
 
        fC = fC + 1;
        
        PCs = interfaceModelAndDiffusionTracker(field,PCs,round(i/fS.FrameSkip)+1);
    end
end