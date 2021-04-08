function outImg = paintOverlapField(WensinkField,pxSize)
%PAINTOVERLAPFIELD Paints ellipses into an input field of zeros such that
%the value in each ellipse corresponds to its index in WensinkField. Used for
%calculating which cells are in contact with each other.

xPxs = ceil(WensinkField.xCells/pxSize);
yPxs = ceil(WensinkField.yCells/pxSize);
majLenPxs = (0.5*WensinkField.aCells*WensinkField.lam)/pxSize;
minLenPxs = (0.5*WensinkField.lam*ones(size(WensinkField.xCells)))/pxSize;
phis = -rad2deg(WensinkField.thetCells);
indices = 1:size(WensinkField.xCells,1);

%These are cludges to account for those very rare occasions when a rod is
%over the edge of the final pixel - just subtract or add one to bring it within
%range.
xPxs(xPxs > round(WensinkField.xWidth/pxSize)) = xPxs(xPxs > round(WensinkField.xWidth/pxSize)) - 1;
yPxs(yPxs > round(WensinkField.yHeight/pxSize)) = xPxs(yPxs > round(WensinkField.yHeight/pxSize)) - 1;
xPxs(xPxs == 0) = 1;
yPxs(yPxs == 0) = 1;

halfWindow = round(max(majLenPxs)) + 3;
fullWindow = 2*halfWindow + 1; %Size of entire cutout

%Prepare image by padding it out.
outImg = zeros(round(WensinkField.xWidth/pxSize),round(WensinkField.yHeight/pxSize));
outImg = [zeros(halfWindow,size(outImg,2));outImg;zeros(halfWindow,size(outImg,2))];
outImg = [zeros(size(outImg,1),halfWindow),outImg,zeros(size(outImg,1),halfWindow)];

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
    if and(and(xPx > -halfWindow,yPx > -halfWindow),and(xPx < (WensinkField.xWidth/pxSize) + halfWindow,yPx < (WensinkField.yHeight/pxSize) + halfWindow))
        minXglob = max(minXbig,1);
        maxXglob = min(maxXbig,floor(WensinkField.xWidth/pxSize) + fullWindow - 1);
        minYglob = max(minYbig,1);
        maxYglob = min(maxYbig,floor(WensinkField.yHeight/pxSize) + fullWindow - 1);
        
        minXloc = minXglob-minXbig+1;
        maxXloc = size(ellipseImg,2) - (maxXbig-maxXglob);
        minYloc = minYglob-minYbig+1;
        maxYloc = size(ellipseImg,1) - (maxYbig-maxYglob);
        
        outImg(minYglob:maxYglob,minXglob:maxXglob) =...
            outImg(minYglob:maxYglob,minXglob:maxXglob)...
            + ellipseImg(minYloc:maxYloc,minXloc:maxXloc)*index;
    end
end

outImg(:,end-halfWindow+1:end) = [];
outImg(:,1:halfWindow) = [];
outImg(end-halfWindow+1:end,:) = [];
outImg(1:halfWindow,:) = [];