function [pD,nD,pDorients,nDorients] = analyseDefects(dx,dy,minimumDist,plotting)
dx4=dx;
dy4=dy;

[nx,ny] = size(dx);

%%%START OF PLOTTING PARAMETERS %%%
nMax = max(nx,ny);
defectScale = nMax/100;
defectThickness = 2;
%%%END OF PLOTTING PARAMETERS%%%

%Find the locations of possible defects (based on sharp edges in the
%orientation field)
se = strel('disk',2,8);
possSpots = imdilate(edge(dy),se);

%Remove objects with holes (defects will be found at the end of line-like
%objects, not ring-like objects)
% CC = bwconncomp(possSpots);
% S = regionprops(CC,'EulerNumber');
% for i = 1:size(S,1)
%     if S(i).EulerNumber < 1
%         possSpots(CC.PixelIdxList{i}) = 0;
%     end
% end

wn=calcs(dx,dy,possSpots);
[pD,nD] = chargearray(nx,ny,wn);

[pD,nD] = annihilateDefects(pD,nD,minimumDist);

[pDorients,nDorients] = findDefectOrients(pD,nD,dx4,dy4,nx,ny,1);

if plotting
    for i = 1:size(nD,1)
        if ~isnan(nDorients(i))
            line([nD(i,2),nD(i,2)+(cosd(nDorients(i)+60)*defectScale*1.6)],[nD(i,1),nD(i,1)+(sind(nDorients(i)+60)*defectScale*1.6)],'Color',[0,0,0],'LineWidth',defectScale/defectThickness);
            line([nD(i,2),nD(i,2)+(cosd(nDorients(i)+60)*defectScale*1.5)],[nD(i,1),nD(i,1)+(sind(nDorients(i)+60)*defectScale*1.5)],'Color',[0,1,1],'LineWidth',defectScale/(defectThickness*2));
            line([nD(i,2),nD(i,2)+(cosd(nDorients(i)+180)*defectScale*1.6)],[nD(i,1),nD(i,1)+(sind(nDorients(i)+180)*defectScale*1.6)],'Color',[0,0,0],'LineWidth',defectScale/defectThickness);
            line([nD(i,2),nD(i,2)+(cosd(nDorients(i)+180)*defectScale*1.5)],[nD(i,1),nD(i,1)+(sind(nDorients(i)+180)*defectScale*1.5)],'Color',[0,1,1],'LineWidth',defectScale/(defectThickness*2));
            line([nD(i,2),nD(i,2)+(cosd(nDorients(i)-60)*defectScale*1.6)],[nD(i,1),nD(i,1)+(sind(nDorients(i)-60)*defectScale*1.6)],'Color',[0,0,0],'LineWidth',defectScale/defectThickness);
            line([nD(i,2),nD(i,2)+(cosd(nDorients(i)-60)*defectScale*1.5)],[nD(i,1),nD(i,1)+(sind(nDorients(i)-60)*defectScale*1.5)],'Color',[0,1,1],'LineWidth',defectScale/(defectThickness*2));
        end
        plotTriangle(nD(i,2),nD(i,1),nDorients(i) + 60,defectScale/1.5,[0,0,0]);
        plotTriangle(nD(i,2),nD(i,1),nDorients(i) + 60,defectScale/2.25,[0,0,1]);
    end
    
    for i = 1:size(pD,1)
        if ~isnan(pDorients(i))
            plotarrow(pD(i,2),pD(i,1),cosd(pDorients(i)),sind(pDorients(i)),[1,0.5,0],defectScale*3,defectScale/defectThickness,defectScale/(defectThickness*2));
        end
        plot(pD(i,2),pD(i,1),'ko','MarkerSize',defectScale/1.5,'MarkerFaceColor','k');
        plot(pD(i,2),pD(i,1),'ro','MarkerSize',defectScale/2.25,'MarkerFaceColor','r');
    end
end