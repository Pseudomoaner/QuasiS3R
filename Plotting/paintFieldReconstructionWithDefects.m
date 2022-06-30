function [] = paintFieldReconstructionWithDefects(trackableData,dS,eF,pDT,fromDefMappings,toDefMappings)

Upsample = 8; %Extent to which the 'design' image should be interpolated to create smoother graphics.
Downsample = 2; %Extent to which the resulting image should be downsized to make the final image size
majorBoost = 0.2; %Extra factor by which to make cells longer for rendering purposes
minorBoost = 0.1; %Extra factor by which to make cells wider for rendering purposes           

%Extract slice representation of defect tracks
defSlice.x = extractTimePoints(fromDefMappings,toDefMappings,pDT,'x');
defSlice.y = extractTimePoints(fromDefMappings,toDefMappings,pDT,'y');
defSlice.phi = extractTimePoints(fromDefMappings,toDefMappings,pDT,'phi');
defSlice.sign = extractTimePoints(fromDefMappings,toDefMappings,pDT,'sparefeat1');

for t = 1:length(trackableData.Length)
    switch dS.colourCells
        case 'Orientation'
            map = colormap('hsv');
            colSet = zeros(size(trackableData.Length{t},1),3);
            for k = 1:size(trackableData.Orientation{t},1)
                bin = ceil((mod(trackableData.Orientation{t}(k)+pi,pi)/pi) * size(map,1));
                colSet(k,:) = map(bin,:);
            end
        case 'Position'
            colSet = zeros(size(trackableData.Length{t},1),3);
            for k = 1:size(trackableData.Length{t},1)
                xFac = trackableData.Centroid{t}(k)/eF.xWidth;
                yFac = trackableData.Centroid{t}(k)/eF.yHeight;
                colSet(k,:) = [xFac,1-(xFac+yFac)/2,yFac];
            end
        case 'Hits' %Choice for debugging purposes - won't necessarily always work with all parameter settings. Assumes a firing and non-firing population.
            maxHit = 5;
            colSet = zeros(size(trackableData.Length{t},1),3);
            for k = 1:size(trackableData.Length{t},1)
                currHits = min(trackableData.Hit{t}(k),maxHit);
                if currHits > 0 %Dead cells are set as magenta shading to gray
                    colFac = (currHits)/(maxHit+1.5);
                    colSet(k,:) = 1-[1,0.5+colFac*0.5,0.05+colFac*0.95];%Orange to white: 1-[1,0.5+colFac*0.5,0.05+colFac*0.95]; Yellow to black: -[1-colFac,0.8-colFac*0.8,0.05-colFac*0.05]; Magenta to black: -[0.9-0.7*colFac,0.2*colFac,0.5-colFac*0.3];%-[0.9,colFac*0.9,0.4+colFac*0.5];
                elseif trackableData.FireRate{t}(k) > 0 %Firing cells are set as light blue
                    colSet(k,:) = 1-([32,128,196]/255);
                elseif trackableData.FireRate{t}(k) == 0 %Unhit cells are set as orange
                    colSet(k,:) = 1-([255, 120, 10]/255);
                end
            end
        case 'None' %Will try to split by population labels
            pop1 = trackableData.Population{t} == 's';
            pop2 = trackableData.Population{t} == 't';
            
            colSet = [0.8*pop1,0.5*pop1 + 0.3*pop2, 0.9*pop2];
    end
    
    backImg = ones(round(eF.yHeight*Upsample),round(eF.xWidth*Upsample));
    xs = trackableData.Centroid{t}(:,1);
    ys = trackableData.Centroid{t}(:,2);
    phis = rad2deg(trackableData.Orientation{t});
    majors = (cos(trackableData.Tilt{t}).*(trackableData.Length{t}-eF.lam) + eF.lam)/2 + majorBoost;
    minors = (eF.lam*ones(size(colSet,1),1)/2)+minorBoost;
    
    %Do separate paintjobs for each of the three colour
    %channels
    imgr = paintEllipse(backImg,xs,ys,majors,minors,phis,colSet(:,1),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    imgg = paintEllipse(backImg,xs,ys,majors,minors,phis,colSet(:,2),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    imgb = paintEllipse(backImg,xs,ys,majors,minors,phis,colSet(:,3),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    
    %Do paintjobs for defects
    phDef = defSlice.sign{t} == 1;
    nhDef = defSlice.sign{t} == -1;
    
    %+1/2
    xs = defSlice.x{t}(phDef);
    ys = defSlice.y{t}(phDef);
    
    imgr = paintCircles(imgr,xs,ys,ones(size(xs))*2,ones(size(xs)),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    imgg = paintCircles(imgg,xs,ys,ones(size(xs))*2,zeros(size(xs)),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    imgb = paintCircles(imgb,xs,ys,ones(size(xs))*2,zeros(size(xs)),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    
    %-1/2
    xs = defSlice.x{t}(nhDef);
    ys = defSlice.y{t}(nhDef);
    phis = defSlice.phi{t}(nhDef);
    
    imgr = paintTriangles(imgr,xs,ys,phis,ones(size(xs))*4,zeros(size(xs)),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    imgg = paintTriangles(imgg,xs,ys,phis,ones(size(xs))*4,zeros(size(xs)),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);
    imgb = paintTriangles(imgb,xs,ys,phis,ones(size(xs))*4,ones(size(xs)),1/Upsample,eF.boundConds,eF.xWidth,eF.yHeight);

    %Smooth to make nicer looking
    imgr = imgaussfilt(imgr,Upsample*eF.lam/4);
    imgg = imgaussfilt(imgg,Upsample*eF.lam/4);
    imgb = imgaussfilt(imgb,Upsample*eF.lam/4);
    
    outImg = cat(3,imgr,imgg,imgb);
    outImg = imresize(outImg,1/Downsample);
    
    Imgpath = sprintf(dS.ImgPath,t-1);
    ImgpathTemp = [dS.imagedirectory, filesep, Imgpath];
    
    imwrite(outImg,ImgpathTemp,'TIFF')
end