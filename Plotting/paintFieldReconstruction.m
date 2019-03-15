function [] = paintFieldReconstruction(trackableData,dS,fS)

resolution = 0.2; %Equivalent to pixel size - how many cell widths each pixel should correspond to

for i = 1:length(trackableData.Length)
    cellImg = zeros(round(fS.fieldHeight/resolution),round(fS.fieldWidth/resolution));
    xs = trackableData.Centroid{i}(:,1);
    ys = trackableData.Centroid{i}(:,2);
    phis = rad2deg(trackableData.Orientation{i});
    majors = (cos(trackableData.Tilt{i}).*(trackableData.Length{i}-fS.lam) + fS.lam)/2;
    minors = repmat(fS.lam/2,size(xs,1),1);
    cellImg = paintEllipse(cellImg,xs,ys,majors,minors,phis,resolution);
    
    Imgpath = sprintf(dS.ImgPath3,i);
    ImgpathTemp = [dS.imagedirectory, filesep, Imgpath];
    
    imwrite(cellImg,ImgpathTemp,'TIFF')
end