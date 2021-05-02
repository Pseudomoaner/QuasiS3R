function expandedProfile = makeExpandedRodProfile(WensinkField,cInd,CDIrange,pxSize)
%Makes an elliptical image of the specified rod at the specified time,
%expanded so you can see which neighbours are within range

%Ellipse parameters

xPx = round(WensinkField.xCells(cInd)/pxSize);
yPx = round(WensinkField.yCells(cInd)/pxSize);
majLenPx = ((WensinkField.aCells(cInd)*WensinkField.lam)+(CDIrange*2))/(2*pxSize);
minLenPx = (WensinkField.lam+(CDIrange*2))/(2*pxSize); %2x actual rod width
phi = rad2deg(WensinkField.thetCells(cInd));

%I will use a slightly more elegant strategy when rewriting this function -
%treat the ellipse as being at the centre of the image, then use circshift
%to move it into the correct location. This will account for the periodic
%boundries neatly.

halfWindow = round(majLenPx) + 3;
fullWindow = 2*halfWindow + 1; %Size of entire cutout (in pixels)

imgXpx = round(WensinkField.xWidth/pxSize);
imgYpx = round(WensinkField.yHeight/pxSize);
expandedProfile = zeros(imgXpx,imgYpx);

%Generate a cut out bit of the coordinate grid for calculating the ellipse
%over.
[xGrid,yGrid] = meshgrid(-halfWindow:halfWindow,-halfWindow:halfWindow);

%Transform x and y grids into canonical coordinates.
xCan = xGrid*cosd(-phi) + yGrid*sind(-phi);
yCan = -xGrid*sind(-phi) + yGrid*cosd(-phi);

geometry = ((xCan.^2)/(majLenPx^2)) + ((yCan.^2)/(minLenPx^2));
ellipseImg = geometry < 1.0; %Image of ellipse in 'small' coordinates

%Draw into middle of image
expandedProfile(round(imgXpx/2) - halfWindow:round(imgXpx/2) + halfWindow,round(imgYpx/2) - halfWindow:round(imgYpx/2) + halfWindow) = ellipseImg;

%Apply circshifts to move ellipse into correct position
expandedProfile = circshift(expandedProfile,yPx-round(imgYpx/2),1);
expandedProfile = circshift(expandedProfile,xPx-round(imgXpx/2),2);

expandedProfile = logical(expandedProfile);