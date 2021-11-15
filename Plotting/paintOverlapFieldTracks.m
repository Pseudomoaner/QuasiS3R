function outImgs = paintOverlapFieldTracks(trackableData,pxSize,width,height)
%Paints ellipses into an input field of zeros such that the value in each
%ellipse corresponds to its track index. Used for calculating which cells
%are in contact with each other.

outImgs = zeros(round(width/pxSize),round(height/pxSize),size(trackableData.Centroid,2));

for i = 1:size(trackableData.Centroid,2)
    xPxs = ceil(trackableData.Centroid{i}(:,1)/pxSize);
    yPxs = ceil(trackableData.Centroid{i}(:,2)/pxSize);
    majLenPxs = (0.5*trackableData.Length{i})/pxSize;
    minLenPxs = (0.5*ones(size(trackableData.Length{i})))/pxSize;
    phis = rad2deg(trackableData.Orientation{i});
    indices = trackableData.TrackIndex{i};
    
    %These are cludges to account for those very rare occasions when a rod is
    %over the edge of the final pixel - just subtract or add one to bring it within
    %range.
    xPxs(xPxs > size(outImgs,2)) = xPxs(xPxs > size(outImgs,2)) - 1;
    yPxs(yPxs > size(outImgs,1)) = xPxs(yPxs > size(outImgs,1)) - 1;
    xPxs(xPxs == 0) = 1;
    yPxs(yPxs == 0) = 1;
    
    halfWindow = round(max(majLenPxs)) + 3;
    fullWindow = 2*halfWindow + 1; %Size of entire cutout
    
    %Prepare image by padding it out.
    thisImg = zeros(round(width/pxSize),round(height/pxSize));
    thisImg = [zeros(halfWindow,size(thisImg,2));thisImg;zeros(halfWindow,size(thisImg,2))];
    thisImg = [zeros(size(thisImg,1),halfWindow),thisImg,zeros(size(thisImg,1),halfWindow)];
    
    for j = 1:size(xPxs,1)
        xPx = xPxs(j);
        yPx = yPxs(j);
        majLenPx = majLenPxs(j);
        minLenPx = minLenPxs(j);
        phi = phis(j);
        index = indices(j);
        
        %Generate a cut out bit of the coordinate grid for calculating the ellipse
        %over.
        [xGrid,yGrid] = meshgrid(xPx-halfWindow:xPx+halfWindow,yPx-halfWindow:yPx+halfWindow);
        
        %Transform x and y grids into canonical coordinates.
        xCan = (xGrid-xPx)*cosd(-phi) + (yGrid-yPx)*sind(-phi);
        yCan = -(xGrid-xPx)*sind(-phi) + (yGrid-yPx)*cosd(-phi);
        
        geometry = ((xCan.^2)/(majLenPx^2)) + ((yCan.^2)/(minLenPx^2));
        ellipseImg = geometry < 1.0; %Image of ellipse in 'small' coordinates
        
        %Do initial drawing
        minXbig = xPx;
        maxXbig = xPx + fullWindow - 1;
        minYbig = yPx;
        maxYbig = yPx + fullWindow - 1;
        
        %Crop the ellipse if it falls outside the image boundary so it fits
        %properly into domain
        if and(and(xPx > -halfWindow,yPx > -halfWindow),and(xPx < (width/pxSize) + halfWindow,yPx < (height/pxSize) + halfWindow))
            minXglob = max(minXbig,1);
            maxXglob = min(maxXbig,floor(width/pxSize) + fullWindow - 1);
            minYglob = max(minYbig,1);
            maxYglob = min(maxYbig,floor(height/pxSize) + fullWindow - 1);
            
            minXloc = minXglob-minXbig+1;
            maxXloc = size(ellipseImg,2) - (maxXbig-maxXglob);
            minYloc = minYglob-minYbig+1;
            maxYloc = size(ellipseImg,1) - (maxYbig-maxYglob);
            
            thisImg(minYglob:maxYglob,minXglob:maxXglob) =...
                thisImg(minYglob:maxYglob,minXglob:maxXglob)...
                + ellipseImg(minYloc:maxYloc,minXloc:maxXloc)*index;
        end
    end
    
    thisImg(:,end-halfWindow+1:end) = [];
    thisImg(:,1:halfWindow) = [];
    thisImg(end-halfWindow+1:end,:) = [];
    thisImg(1:halfWindow,:) = [];
    
    outImgs(:,:,i) = thisImg;
end