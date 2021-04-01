function pHalfOri =GetphalfDefectOrientationNew(xphdcentre,yphdcentre,NemAngleListCW,XPositionNemList,YPositionNemList)

XPosCorVectList=xphdcentre-XPositionNemList; %position to core vector, x component
YPosCorVectList=yphdcentre-YPositionNemList; %position to core vector, y component

phiList = atan2(YPosCorVectList,XPosCorVectList);
FT = sum(exp(1i*phiList).*cos(2*NemAngleListCW)); %Fourier transform of the positional and orientational data
pHalfOri = rad2deg(-angle(FT));