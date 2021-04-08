function [negDefCents,negDefOris,posDefCents,posDefOris] = measureDefects(trackableData,fieldSettings,procSettings)

negDefCents = cell(size(trackableData.Centroid,1),1);
posDefCents = cell(size(trackableData.Centroid,1),1);
posDefOris = cell(size(trackableData.Centroid,1),1);
negDefOris = cell(size(trackableData.Centroid,1),1);

plotting = false; %Set to true for debugging
imgDir = 'C:\Users\olijm\Desktop\TmpImgs';

for i = 1:size(trackableData.Centroid,2)
    imgDat = fieldReconst(trackableData,fieldSettings,procSettings,i);
    oriDat = findImageOrients(imgDat,procSettings.tensorSize/procSettings.pixSize);
    
    oriX = cos(oriDat);
    oriY = -sin(oriDat);
    
    if plotting
        figure(1)
        cla
        imshow(imgDat,[])
        hold on
    end
    
    [posDefCents{i},negDefCents{i},posDefOris{i},negDefOris{i}] = analyseDefects(oriX,oriY,0,plotting);
    
    if plotting
        [nr,nc] = size(imgDat);
        axis([1,nr,1,nc])
        export_fig(fullfile(imgDir,sprintf('Frame_%04d.tif',i)),'-nocrop')
    end
    
    posDefCents{i} = flip(posDefCents{i},2)*procSettings.pixSize;
    negDefCents{i} = flip(negDefCents{i},2)*procSettings.pixSize;
end