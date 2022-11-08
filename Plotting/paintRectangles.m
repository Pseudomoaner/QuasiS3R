function outImg = paintRectangles(inImg,xs,ys,thetas,edgeLens,edgeWids,intensities,pxSize,boundaryConds,width,height)
%PAINTTRIANGLES paints a set of triangles of the given position, scale, 
%colour and orientation.

xPxs = ceil(xs/pxSize);
yPxs = ceil(ys/pxSize);

%These are cludges to account for those very rare occasions when a defect is
%over the edge of the final pixel - just subtract or add one to bring it within
%range.
xPxs(xPxs > size(inImg,2)) = xPxs(xPxs > size(inImg,2)) - 1;
yPxs(yPxs > size(inImg,1)) = xPxs(yPxs > size(inImg,1)) - 1;
xPxs(xPxs == 0) = 1;
yPxs(yPxs == 0) = 1;

offsets = (edgeLens+edgeWids)/2;
offsetPxs = ceil(offsets/pxSize);

halfWindow = round(max(offsetPxs)*2) + 3;
fullWindow = 2*halfWindow + 1; %Size of entire cutout

%Prepare image by padding it out.
outImg = inImg;
if strcmp(boundaryConds,'periodic')
    outImg = [outImg(end-halfWindow+1:end,:);outImg;outImg(1:halfWindow,:)];
    outImg = [outImg(:,end-halfWindow+1:end),outImg,outImg(:,1:halfWindow)];
else %This isn't really super necessary, but means you have to do less coding later as it makes periodic and non-periodic image coordinate spaces the same.
    outImg = [zeros(halfWindow,size(outImg,2));outImg;zeros(halfWindow,size(outImg,2))];
    outImg = [zeros(size(outImg,1),halfWindow),outImg,zeros(size(outImg,1),halfWindow)];
end

for i = 1:size(xs,1)
    xPx = xPxs(i);
    yPx = yPxs(i);
    theta = thetas(i);
    edgeLenPx = edgeLens(i)/pxSize;
    edgeWidPx = edgeWids(i)/pxSize;
    
    %Generate a cut out bit of the coordinate grid for calculating the
    %rectangle over.
    [xGrid,yGrid] = meshgrid(-halfWindow:halfWindow,-halfWindow:+halfWindow);
    
    xRot = cosd(theta) * xGrid - sind(theta) * yGrid;
    yRot = sind(theta) * xGrid + cosd(theta) * yGrid;
    
    rectImg = and(and(xRot < edgeLenPx/2,xRot > -edgeLenPx/2),and(yRot < edgeWidPx/2,yRot > -edgeWidPx/2));
    
    %Do initial drawing
    minXbig = xPx;
    maxXbig = xPx + fullWindow - 1;
    minYbig = yPx;
    maxYbig = yPx + fullWindow - 1;
    
    if strcmp(boundaryConds,'periodic')
        %Don't need to worry about triangles that are outside the actual field,
        %as the periodic boundary conditions will move them back in.
        existingImg = outImg(minYbig:maxYbig,minXbig:maxXbig);
        existingImg(rectImg) = intensities(i);
        outImg(minYbig:maxYbig,minXbig:maxXbig) = existingImg;
    else %Non-periodic boundary conditions. Here you do have to worry about rods outside the main field.
        if and(and(xPx > -halfWindow,yPx > -halfWindow),and(xPx < (width/pxSize) + halfWindow,yPx < (height/pxSize) + halfWindow))
            minXglob = max(minXbig,1);
            maxXglob = min(maxXbig,floor(width/pxSize) + fullWindow - 1);
            minYglob = max(minYbig,1);
            maxYglob = min(maxYbig,floor(height/pxSize) + fullWindow - 1);
            
            minXloc = minXglob-minXbig+1;
            maxXloc = size(rectImg,2) - (maxXbig-maxXglob);
            minYloc = minYglob-minYbig+1;
            maxYloc = size(rectImg,1) - (maxYbig-maxYglob);

            existingImg = outImg(minYglob:maxYglob,minXglob:maxXglob);
            existingImg(rectImg) = intensities(i);
            outImg(minYglob:maxYglob,minXglob:maxXglob) = existingImg;
        end
    end
end

if strcmp(boundaryConds,'periodic')
    %Cut out padded edges and move them around to do wrap around
    %x-dimension first...
    leftChunk1 = outImg(:,halfWindow+1:fullWindow);
    changeInds1 = leftChunk1 < 0.999;
    jointChunk1 = outImg(:,end-halfWindow:end);
    jointChunk1(changeInds1) = leftChunk1(changeInds1);
    outImg(:,halfWindow+1:fullWindow) = jointChunk1;
    
    rightChunk2 = outImg(:,end-fullWindow:end-halfWindow-2);
    changeInds2 = rightChunk2 < 0.999;
    jointChunk2 = outImg(:,1:halfWindow);
    jointChunk2(changeInds2) = rightChunk2(changeInds2);
    outImg(:,end-fullWindow:end-halfWindow-2) = jointChunk2;
    
    %Then y
    bottomChunk3 = outImg(halfWindow+1:fullWindow,:);
    changeInds3 = bottomChunk3 < 0.999;
    jointChunk3 = outImg(end-halfWindow:end,:);
    jointChunk3(changeInds3) = bottomChunk3(changeInds3);
    outImg(halfWindow+1:fullWindow,:) = jointChunk3;
    
    topChunk4 = outImg(end-fullWindow:end-halfWindow-2,:);
    changeInds4 = topChunk4 < 0.999;
    jointChunk4 = outImg(1:halfWindow,:);
    jointChunk4(changeInds4) = topChunk4(changeInds4);
    outImg(end-fullWindow:end-halfWindow-2,:) = jointChunk4;
end

outImg(:,end-halfWindow+1:end) = [];
outImg(:,1:halfWindow) = [];
outImg(end-halfWindow+1:end,:) = [];
outImg(1:halfWindow,:) = [];
