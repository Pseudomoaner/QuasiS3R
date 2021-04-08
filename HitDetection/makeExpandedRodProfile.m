function expandedProfile = makeExpandedRodProfile(WensinkField,cInd,CDIrange,pxSize)
%Makes an elliptical image of the specified rod at the specified time,
%expanded so you can see which neighbours are within range

xPx = round(WensinkField.xCells(cInd)/pxSize);
yPx = round(WensinkField.yCells(cInd)/pxSize);
majLenPx = ((WensinkField.aCells(cInd)*WensinkField.lam)+(CDIrange*2))/(2*pxSize);
minLenPx = (WensinkField.lam+(CDIrange*2))/(2*pxSize); %2x actual rod width
phi = -rad2deg(WensinkField.thetCells(cInd));

halfWindow = round(majLenPx) + 3;
fullWindow = 2*halfWindow + 1; %Size of entire cutout

expandedProfile = zeros(round(WensinkField.xWidth/pxSize),round(WensinkField.yHeight/pxSize));
expandedProfile = [zeros(halfWindow,size(expandedProfile,2));expandedProfile;zeros(halfWindow,size(expandedProfile,2))];
expandedProfile = [zeros(size(expandedProfile,1),halfWindow),expandedProfile,zeros(size(expandedProfile,1),halfWindow)];

%Generate a cut out bit of the coordinate grid for calculating the ellipse
%over.
[xGrid,yGrid] = meshgrid(xPx-halfWindow:xPx+halfWindow,yPx-halfWindow:yPx+halfWindow);

%Transform x and y grids into canonical coordinates.
xCan = (xGrid-xPx)*cosd(-phi) + (yGrid-yPx)*sind(-phi);
yCan = -(xGrid-xPx)*sind(-phi) + (yGrid-yPx)*cosd(-phi);

geometry = ((xCan.^2)/(majLenPx^2)) + ((yCan.^2)/(minLenPx^2));
ellipseImg = geometry < 1.0; %Image of ellipse in 'small' coordinates

%Do drawing
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
    
    expandedProfile(minYglob:maxYglob,minXglob:maxXglob) = ellipseImg(minYloc:maxYloc,minXloc:maxXloc);
end

expandedProfile(:,end-halfWindow+1:end) = [];
expandedProfile(:,1:halfWindow) = [];
expandedProfile(end-halfWindow+1:end,:) = [];
expandedProfile(1:halfWindow,:) = [];

expandedProfile = logical(expandedProfile);