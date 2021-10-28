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

switch patchSettings.patchType
    case 'Circle'
        %Make an initial guess at the radius of the circle by assuming the field is
        %uniformly populated
        testRad = sqrt(patchSettings.fraction*inField.xWidth*inField.yHeight/pi);
        actualFrac = findFracInCircle(inField,testRad);

        tooSmall = actualFrac < patchSettings.fraction - (patchSettings.tol*patchSettings.fraction);
        tooBig = actualFrac > patchSettings.fraction + (patchSettings.tol*patchSettings.fraction);

        stepCnt = 0;

        while tooBig || tooSmall
            %This bit of code terminates the while loop if it's been running for
            %too long with no improvement
            stepCnt = stepCnt + 1;
            if stepCnt > 10
                warning('Maximum number of steps in the patch radius calculation exceeded. Using best approximation...')
                break
            end

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
    case 'Voronoi'
        noSeeds = round(inField.xWidth*inField.yHeight*patchSettings.seedDensity);
        noSeeds1 = round(noSeeds * patchSettings.seedFrac);
        noSeeds2 = noSeeds - noSeeds1;
        seedsX = rand(noSeeds,1)*inField.xWidth;
        seedsY = rand(noSeeds,1)*inField.yHeight;
        tmpSeeds = [zeros(noSeeds2,1);ones(noSeeds1,1)];
        popLabels = tmpSeeds(randperm(noSeeds)); %0 indicates pop 2, 1 pop 1
        
        %Create a tiling of these initial seeds so the resulting Voronoi
        %diagram is periodic
        seedsX = [seedsX;seedsX;seedsX;seedsX - inField.xWidth;seedsX - inField.xWidth;seedsX - inField.xWidth;seedsX + inField.xWidth;seedsX + inField.xWidth;seedsX + inField.xWidth];
        seedsY = [seedsY;seedsY - inField.yHeight;seedsY + inField.yHeight;seedsY;seedsY - inField.yHeight;seedsY + inField.yHeight;seedsY;seedsY - inField.yHeight;seedsY + inField.yHeight];
        popLabels = repmat(popLabels,9,1);

        DT = delaunayTriangulation([seedsX,seedsY]);
        nearestSeeds = nearestNeighbor(DT,inField.xCells,inField.yCells);

        for i = 1:size(inField.xCells,1)
            if popLabels(nearestSeeds(i))
                inField.popCells(i) = patchSettings.popLabel;
                inField.fCells(i) = patchSettings.force;
                inField.cCells(i,:) = patchSettings.colour;
                inField.rCells(i) = patchSettings.reversalRate;
                inField.fireCells(i) = patchSettings.fireRate;
            end
        end

        outField = inField;
end
end

function popFrac = findFracInCircle(inField,rad)
    %Finds the actual fraction of rods inside the circle centred on
    %(inField.xWidth/2,inField.yHeight/2) of radius rad.
    xs = inField.xCells - inField.xWidth/2;
    ys = inField.yCells - inField.yHeight/2;
    
    popSize = sum(sqrt(xs.^2 + ys.^2)<rad);
    popFrac = popSize/size(inField.xCells,1);
end