function outField = makePatch(inField,patchSettings)
%MAKEPATCH creates a localised patch of a particular cell type in the
%centre of the simulation domain.
%
%   INPUTS:
%       -inField: The WensinkField simulation you want to apply the patch
%       to
%       -patchSettings: Structure defining the size and properties of the
%       patch
%
%   OUTPUTS:
%       -outField: The WensinkField simulation with the circular patch
%       applied.
%
%   Author: Oliver J. Meacock, (c) 2021

%Make an initial guess at the radius of the circle by assuming the field is
%uniformly populated
testRad = sqrt(patchSettings.fraction*inField.xWidth*inField.yHeight/pi);
actualFrac = findFracInCircle(inField,testRad);

tooSmall = actualFrac < patchSettings.fraction - (patchSettings.tol*patchSettings.fraction);
tooBig = actualFrac > patchSettings.fraction + (patchSettings.tol*patchSettings.fraction);

while tooBig || tooSmall
    fracScale = patchSettings.fraction/actualFrac;
    
    testRad = testRad*sqrt(fracScale);
    actualFrac = findFracInCircle(inField,testRad);
    
    tooSmall = actualFrac < patchSettings.fraction - (patchSettings.tol*patchSettings.fraction);
    tooBig = actualFrac > patchSettings.fraction + (patchSettings.tol*patchSettings.fraction);
end

xs = inField.xCells - inField.xWidth/2;
ys = inField.yCells - inField.yHeight/2;

changeInds = sqrt(xs.^2 + ys.^2)<testRad;

inField.fCells(changeInds) = patchSettings.force;
inField.cCells(changeInds,:) = repmat(patchSettings.colour,sum(changeInds),1);
inField.rCells(changeInds) = patchSettings.reversalRate;
inField.fireCells(changeInds) = patchSettings.fireRate;
inField.popCells(changeInds) = repmat(patchSettings.popLabel,sum(changeInds),1);

outField = inField;
end

function popFrac = findFracInCircle(inField,rad)
    %Finds the actual fraction of rods inside the circle centred on
    %(inField.xWidth/2,inField.yHeight/2) of radius rad.
    xs = inField.xCells - inField.xWidth/2;
    ys = inField.yCells - inField.yHeight/2;
    
    popSize = sum(sqrt(xs.^2 + ys.^2)<rad);
    popFrac = popSize/size(inField.xCells,1);
end