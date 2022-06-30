function cellImg = fieldReconst(trackableData,endField,pS,frameInd)

cellImg = zeros(round(endField.yHeight/pS.pixSize),round(endField.xWidth/pS.pixSize));
xs = trackableData.Centroid{frameInd}(:,1);
ys = trackableData.Centroid{frameInd}(:,2);
phis = rad2deg(trackableData.Orientation{frameInd});
majors = (cos(trackableData.Tilt{frameInd}).*(trackableData.Length{frameInd}-endField.lam) + endField.lam)/2;
minors = repmat(endField.lam/2,size(xs,1),1);
cellImg = paintEllipse(cellImg,xs,ys,majors,minors,phis,ones(size(xs)),pS.pixSize,endField.boundConds,endField.xWidth,endField.yHeight);