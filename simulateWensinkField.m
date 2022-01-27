function [PCs,field,hitNos,contactSet] = simulateWensinkField(startField,fS,dS,contactFind)

PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;

contactSet = {};
hitNos.new = zeros(ceil(fS.motileSteps/fS.FireSkip),1);
hitNos.tot = zeros(ceil(fS.motileSteps/fS.FireSkip),1);

%Actual simulation
for i = 1:fS.motileSteps
    fprintf('Main frame is %i of %i\n',i,fS.motileSteps)
    field = field.stepModel(fS.motiledt,fS.growthRate,fS.divThresh,fS.postDivMovement,dS.colourCells,'main');
    
    %Simulate firing events
    if rem(i,fS.FireSkip) == 0 && sum(field.fireCells) > 0 %Second condition prevents you from bothering with expensive calculations if no cell can fire
        [field,hitNos.new(round(i/fS.FireSkip)),hitNos.tot(round(i/fS.FireSkip))] = field.calculateHits(fS.FireSkip*fS.motiledt);
        field = field.killCells();
        
        if contactFind
            contactSet{round(i/fS.FireSkip)} = field.calculateContacts();
        end
    end
    
    %Output visualisation
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