function [PCs,field] = simulateWensinkFieldBurnIn(startField,fS,dS)
%Separate simulation to allow simulation to settle into passive configuration without cell motility, growth or tilting taking place.
PCs = [];
field = startField;

PCs = interfaceModelAndDiffusionTracker(field,PCs,1);

fC = 0;
for i = 1:fS.burnInSteps
    fprintf('Frame is %i\n',i)
    
    field = field.stepModel(fS.burnIndt,fS.f0,inf,0,fS.divThresh,fS.postDivMovement,fS.colJigRate,dS.colourCells);
   
    if rem(i,fS.FrameSkip) == 0
%         field.drawField(dS.posVec,dS.colourCells,fS.maxX);
%         pause(0.01)
        if dS.saveFrames
            imPath1 = sprintf(dS.ImgPath1,fC);
            fullImPath = [dS.imagedirectory, filesep, imPath1];
            
            export_fig(fullImPath,'-tif','-m1')
            reopenimage = imread(fullImPath);
            cropped=imcrop(reopenimage,dS.cropRect);
            imwrite(cropped,fullImPath,'tif')
            
            disp([fullImPath,' saved.'])
            
            cla
        end
 
        fC = fC + 1;
        
        PCs = interfaceModelAndDiffusionTracker(field,PCs,round(i/fS.FrameSkip)+1);
    end
end