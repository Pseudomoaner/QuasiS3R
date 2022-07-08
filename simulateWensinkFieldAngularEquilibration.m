function [PCs,field] = simulateWensinkFieldAngularEquilibration(startField,fS,dS)
%Separate simulation to allow simulation to settle into active (rather than previous passive) configuration without cell tilting taking place.
PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;

horizontals = abs(cos(PCs.Orientation{end})) > sqrt(2)/2; %Marks cells that are closer to horizontal than vertical
angRat = max(sum(horizontals),sum(~horizontals))/min(sum(horizontals),sum(~horizontals));

stepInd = 1;

while angRat > 1.25 %Will consider equilibrated once this condition is satisfied
    
    field = field.stepModel(fS.motiledt,fS.growthRate,fS.divThresh,fS.postDivMovement,dS.colourCells,'Initial');
    
    if rem(stepInd,fS.FrameSkip) == 0
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
        
        PCs = interfaceModelAndDiffusionTracker(field,PCs,round(stepInd/fS.FrameSkip)+1);
        
        horizontals = abs(cos(PCs.Orientation{end})) > sqrt(2)/2; %Marks cells that are closer to horizontal than vertical
        angRat = max(sum(horizontals),sum(~horizontals))/min(sum(horizontals),sum(~horizontals));
        
        fprintf('Angular ratio: %f (target 1.25)\n',angRat)
    end
    stepInd = stepInd + 1;
end