classdef WensinkField
    properties
        xWidth %Width of field
        yHeight %Height of field
        zDepth %Depth of field
        
        xCells %x positions of centers of mass of cells in field
        yCells %Likewise for y positions
        zCells %And for z positions
        thetCells %Likewise for orientations of cells in xy-plane
        phiCells %Tilts of all cells (how they are aligned relative to the xy-plane in degrees)
        aCells %Likewise for aspect ratios of cells
        nCells %Likewise for number of segments in cells
        uCells %Likewise for orientation vector of cells
        lCells %Likewise for separation between Yukawa segments of cells.
        fCells %Likewise for the extra force added in the direction of motion (The 'self propelled' bit)
        rCells %Likewise for the reversal rate of each cell (probability of reversal in a single time unit)
        cCells %Likewise for color vector of cells.
        hitCells %Likewise for number of hits each cell has experienced
        fireCells %Likewise for firing rate of each cell
        popCells %Likewise for population label for each cell
        
        xBarr %x position for the barrier rods (non-motile, but still cause steric effects)
        yBarr
        zBarr
        thetBarr
        phiBarr
        aBarr %The approach here will be to model the boundary segments as rods of aspect ratio a = 1. 
        nBarr %This is not the most elegant solution, but allows us to recycle the code for rod-rod interactions more easily.
        lBarr
        
        cellDistsBiophys %Distance between all cells. Indexed in same way as obj.cells.
        cellDistsContact %Distance between all cells. Used for CDI calculations, rather than biophysics.
        distThresh %Distance between cells, beyond which biophysical interactions are too weak to be relevant. 
        contactThresh %Distance between cells, beyond which CDI cannot occur.

        U0 %Potential Amplitude
        f0 %Stoksian friction coefficient
        frictionAnisotropy %How isotropic the friction applied to each rod should be - 0 implies isotropic, 1 implies the default anisotropy (from the original Wensink model)
        zElasticity %Elasticity of the confining material in the z-direction
        lam %Screening length of Yukawa segements
        contactRange %Range of CDI system
        killType %Way dead cells are dealt with (lyse or husk)
        killThresh %How many hits must accumulate for a sensitive cell to be killed
        hitRateType %Way CDI hits are dealt out to neighbours (distributed or constant)

        boundConds %Boundary conditions (periodic or none)
        boundaryDesign %Image containing an image of the design of the confinement
        resUp %Difference in resolution between boundaryDesign image and actual simulation domain
        colJigRate %How quickly colours should 'jiggle' noisily (in HSV space). Values above 0 can be useful for visualising lineages of dividing cells.
        
        gpuAvailable %Whether there is a GPU available in the current system
        compiled %Whether there is a compiled .mex version of the potential gradient functions available
    end
    methods
        function obj = WensinkField(fS)
            inputCheck = isfield(fS,{'xWidth','yHeight','zDepth','U0','lam','boundaryConditions','f0','zElasticity','contactRange','killType','hitRateType','colJigRate','killThresh','frictionAnisotropy'});

            %Check and store domain size settings
            if ~inputCheck(1) || ~inputCheck(2) %Must have the domain size
                error('fieldSettings must include xWidth and yHeight as fields')
            else
                validateattributes(fS.xWidth,{'numeric'},{'scalar','positive'})
                validateattributes(fS.yHeight,{'numeric'},{'scalar','positive'})
                obj.xWidth = fS.xWidth;
                obj.yHeight = fS.yHeight;

                if inputCheck(3)
                    validateattributes(fS.zDepth,{'numeric'},{'scalar','positive'})
                    obj.zDepth = fS.zDepth;
                else
                    obj.zDepth = 20; %Default value
                end
            end            
            
            %Check U0
            if inputCheck(4)
                validateattributes(fS.U0,{'numeric'},{'scalar','positive'})
                obj.U0 = fS.U0;
            else
                obj.U0 = 250;
            end
            
            %Check lambda
            if inputCheck(5)
                validateattributes(fS.lam,{'numeric'},{'scalar','positive'})
                obj.lam = fS.lam;
            else
                obj.lam = 1;
            end

            %Check boundary conditions
            if inputCheck(6)
                validateattributes(fS.boundaryConditions,{'char'},{'scalartext'})
                if ~any(strcmp({'none','periodic'},fS.boundaryConditions))
                    error('Expected boundary conditions to be specified as either "none" or "periodic".')
                else
                    obj.boundConds = fS.boundaryConditions;
                end
            else
                obj.boundConds = 'periodic';
            end

            %Check friction coefficient
            if inputCheck(7)
                validateattributes(fS.f0,{'numeric'},{'scalar','positive'})
                obj.f0 = fS.f0;
            else
                obj.f0 = 1;
            end

            %Check z-dimensional elasticity
            if inputCheck(8)
                validateattributes(fS.zElasticity,{'numeric'},{'scalar','positive'})
                obj.zElasticity = fS.zElasticity;
            else
                obj.zElasticity = inf;
            end

            %Check firing range
            if inputCheck(9)
                validateattributes(fS.contactRange,{'numeric'},{'scalar','positive'})
                obj.contactRange = fS.contactRange;
            else
                obj.contactRange = 2;
            end

            %Check killing type
            if inputCheck(10)
                validateattributes(fS.killType,{'char'},{'scalartext'})
                if ~any(strcmp({'lyse','husk'},fS.killType))
                    error('Expected killing type to be specified as either "lyse" or "husk".')
                else
                    obj.killType = fS.killType;
                end
            else
                obj.killType = 'husk';
            end

            %Check killing type
            if inputCheck(11)
                validateattributes(fS.hitRateType,{'char'},{'scalartext'})
                if ~any(strcmp({'distributed','constant'},fS.hitRateType))
                    error('Expected hit rate type to be specified as either "distributed" or "constant".')
                else
                    obj.hitRateType = fS.hitRateType;
                end
            else
                obj.hitRateType = 'distributed';
            end

            %Check colour drift rate
            if inputCheck(12)
                validateattributes(fS.colJigRate,{'numeric'},{'scalar','positive'})
                obj.colJigRate = fS.colJigRate;
            else
                obj.colJigRate = 0;
            end

            %Check killing threshold
            if inputCheck(13)
                validateattributes(fS.killThresh,{'numeric'},{'scalar','positive','integer'})
                obj.killThresh = fS.killThresh;
            else
                obj.killThresh = 1;
            end

            if inputCheck(14)
                validateattributes(fS.frictionAnisotropy,{'numeric'},{'scalartext','positive'})
                obj.frictionAnisotropy = fS.frictionAnisotropy;
            else
                obj.frictionAnisotropy = 1;
            end

            %Check to see if GPU and/or .mex files are available for
            %calculating the potential between rods
            try
                gpuArray(1);
                obj.gpuAvailable = true;
            catch
                obj.gpuAvailable = false;
            end
            
            locPath = mfilename('fullpath');
            locPath = locPath(1:end-13); %Path without the uneccessary '\WensinkField' on the end
            
            if exist(fullfile(locPath,'PotentialCalculations',['mexCalcEnergyGradientsPeriodic.',mexext]),'file') && strcmp(obj.boundConds,'periodic')
                obj.compiled = true;
            elseif exist(fullfile(locPath,'PotentialCalculations',['mexCalcEnergyGradients.',mexext]),'file') && strcmp(obj.boundConds,'none')
                obj.compiled = true;
            else
                obj.compiled = false;
            end
        end
        
        function obj = populateField(obj,barrierSettings,cellSettings,areaFrac)
            %Populates the field with a number of randomly positioned cells of the given force generation and a number of static 'barrier' rods.
            
            switch barrierSettings.type
                case 'none'
                    obj.xBarr = [];
                    obj.yBarr = [];
                    obj.zBarr = [];
                    obj.thetBarr = [];
                    obj.phiBarr = [];
                    obj.aBarr = [];
                    obj.nBarr = [];
                    obj.lBarr = [];
                    obj.boundaryDesign = zeros(round(obj.yHeight),round(obj.xWidth));
                    obj.resUp = 1;
                case 'loaded'
                    obj.xBarr = barrierSettings.CPs(:,1);
                    obj.yBarr = barrierSettings.CPs(:,2);
                    obj.zBarr = ones(size(obj.xBarr)) * obj.zDepth/2;
                    obj.thetBarr = zeros(size(obj.xBarr));
                    obj.phiBarr = zeros(size(obj.xBarr));
                    obj.aBarr = ones(size(obj.xBarr));
                    obj.nBarr = ones(size(obj.xBarr));
                    obj.lBarr = zeros(size(obj.xBarr));
                    obj.boundaryDesign = barrierSettings.fieldImg;
                    obj.resUp = barrierSettings.resUp;
            end     
            
            switch cellSettings.type
                case 'singleCell'
                    obj.xCells = obj.xWidth/2;
                    obj.yCells = obj.yHeight/2;
                    obj.zCells = obj.zDepth/2;
                    obj.thetCells = pi/2;
                    obj.phiCells = 0;
                    obj.aCells = cellSettings.a;
                    obj.fCells = cellSettings.f;
                    obj.rCells = cellSettings.r;
                    obj.cCells = [1,0,0];
                    obj.hitCells = 0;
                    obj.popCells = 's';
                    obj.fireCells = cellSettings.fire;
                case 'doubleCell'
                    obj.xCells = [obj.xWidth/2+1;obj.xWidth/2-cellSettings.a2];
                    obj.yCells = [obj.yHeight/2;obj.yHeight/2];
                    obj.zCells = [obj.zDepth/2;obj.zDepth/2];
                    obj.thetCells = [pi/2;0];
                    obj.phiCells = [0;0];
                    obj.aCells = [cellSettings.a1;cellSettings.a2];
                    obj.fCells = [cellSettings.f1;cellSettings.f2];
                    obj.rCells = [cellSettings.r1;cellSettings.r2];
                    obj.cCells = [1,0,0;0,1,0];
                    obj.hitCells = [0;0];
                    obj.popCells = [cellSettings.pop1;cellSettings.pop2];
                    obj.fireCells = [cellSettings.fire1;cellSettings.fire2];
                case 'LatticedXYCells' %Start with an initially ordered lattice of cells (with random up/down orientations)
                    obj = obj.overlayLattice(areaFrac,cellSettings.a);
                    obj.zCells = ones(size(obj.xCells)) * obj.zDepth/2;
                    obj.phiCells = abs(randn(size(obj.xCells)))*0.001;
                    obj.aCells = cellSettings.a * ones(size(obj.xCells));
                    obj.fCells = cellSettings.f * ones(size(obj.xCells));
                    obj.rCells = cellSettings.r * ones(size(obj.xCells));
                    obj.cCells = repmat(cellSettings.c,size(obj.xCells,1),1);
                    obj.hitCells = zeros(size(obj.xCells));
                    obj.popCells = repmat(cellSettings.pop,size(obj.xCells,1),1);
                    obj.fireCells = cellSettings.fire * ones(size(obj.xCells));
                case 'LatticedXYCellsNoBarrier' %Uses a less buggy method for inserting cells into the grid - useful if you don't need to worry about barriers
                    totArea = (obj.xWidth * obj.yHeight) - (sum(obj.boundaryDesign(:))/(obj.resUp^2));
                    singleArea = (obj.lam^2 * (cellSettings.a - 1)) + (pi * (obj.lam/2)^2);
                    
                    tgtRodNo = round(totArea*areaFrac/singleArea); %Total number of rods you want inserted into simulation domain
                    
                    %Space rods ~ 1 unit apart along their long axes
                    noX = floor(obj.xWidth / (cellSettings.a + 0.5)); %Total number of cells laid end-to-end
                    xSpace = obj.xWidth/noX;
                    
                    %And infer the necesssary y-spacing from these
                    %calculations
                    noY = ceil(tgtRodNo/noX);
                    ySpace = obj.yHeight/noY;
                    
                    %Place rods into the grid until you reach the target
                    %number (will leave a small patch in bottom right of
                    %domain without rods - will disappear during approach
                    %to steady-state)
                    currNo = 1;
                    for i = 1:noX
                        for j = 1:noY
                            obj.xCells(currNo,1) = xSpace*i - (xSpace/2);
                            obj.yCells(currNo,1) = ySpace*j - (ySpace/2);
                            if currNo == tgtRodNo
                                break
                            end
                            currNo = currNo + 1;
                        end
                    end
                    obj.xCells = obj.xCells + rand(tgtRodNo,1)*0.0001;
                    obj.yCells = obj.yCells + rand(tgtRodNo,1)*0.0001;
                    obj.zCells = ones(tgtRodNo,1) * obj.zDepth/2;
                    obj.thetCells = (rand(tgtRodNo,1)>0.5)*pi;
                    obj.phiCells = abs(randn(tgtRodNo,1))*0.001;
                    obj.aCells = cellSettings.a * ones(tgtRodNo,1);
                    obj.fCells = cellSettings.f * ones(tgtRodNo,1);
                    obj.rCells = cellSettings.r * ones(tgtRodNo,1);
                    obj.cCells = rand(tgtRodNo,3);
                    obj.hitCells = zeros(tgtRodNo,1);
                    obj.popCells = repmat(cellSettings.pop,size(obj.xCells,1),1);
                    obj.fireCells = cellSettings.fire * ones(tgtRodNo,1);
                case 'LatticedXYCellsTwoPops' %Start with an initially ordered lattice of cells (with random up/down orientations)
                    obj = obj.overlayLattice(areaFrac,cellSettings.a1); %Because of the way overlayLattice works, ensure that a1 is larger than or equal to a2.
                    obj.zCells = ones(size(obj.xCells)) * obj.zDepth/2;
                    obj.phiCells = abs(randn(size(obj.xCells)))*0.001;
                    
                    %Create randomised initial population mix
                    type = zeros(size(obj.xCells,1),1);
                    type(1:round(cellSettings.popFrac*size(obj.xCells,1))) = 1;
                    type = logical(type(randperm(size(obj.xCells,1))));
                    obj.aCells = zeros(size(obj.xCells,1),1);
                    obj.fCells = zeros(size(obj.xCells,1),1);
                    obj.rCells = zeros(size(obj.xCells,1),1);
                    obj.cCells = zeros(size(obj.xCells,1),3);
                    obj.hitCells = zeros(size(obj.xCells));
                    obj.fireCells = zeros(size(obj.xCells));
                    obj.popCells = repmat(' ',size(obj.xCells,1),1);
                    obj.aCells(type) = cellSettings.a1 * ones(sum(type),1);
                    obj.fCells(type) = cellSettings.f1 * ones(sum(type),1);
                    obj.rCells(type) = cellSettings.r1 * ones(sum(type),1);
                    obj.cCells(type,:) = repmat(cellSettings.c1,sum(type),1);
                    obj.popCells(type) = cellSettings.pop1;
                    obj.fireCells(type) = cellSettings.fire1 * ones(sum(type),1);
                    obj.aCells(~type) = cellSettings.a2 * ones(size(obj.xCells,1) - sum(type),1);
                    obj.fCells(~type) = cellSettings.f2 * ones(size(obj.xCells,1) - sum(type),1);
                    obj.rCells(~type) = cellSettings.r2 * ones(size(obj.xCells,1) - sum(type),1);
                    obj.cCells(~type,:) = repmat(cellSettings.c2,sum(~type),1);
                    obj.popCells(~type) = cellSettings.pop2;
                    obj.fireCells(~type) = cellSettings.fire2 * ones(sum(~type),1);
            end
            
            [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
            obj.uCells = [cos(obj.thetCells).*cos(obj.phiCells),sin(obj.thetCells).*cos(obj.phiCells),sin(obj.phiCells)];
        end
        
        function obj = overlayLattice(obj,areaFrac,aspRat)
            %Overlays a lattice of rods on top of the region defined by the
            %design image, and adjusts density until it matches the desired
            %area fraction. Assumes all rods are identical in size.
            totArea = (obj.xWidth * obj.yHeight) - (sum(obj.boundaryDesign(:))/(obj.resUp^2));
            singleArea = (obj.lam^2 * (aspRat - 1)) + (pi * (obj.lam/2)^2);
            
            tgtRodNo = round(totArea*areaFrac/singleArea);
            
            %Approach will be to do a lattice that is slightly too loose,
            %and then gradually tighten it until the target no. of rods fit within
            %the bounded domain.
            tryYSpace = sqrt(singleArea/(areaFrac*(aspRat+1))); %Spacing of lattice in y-direction
            tryXSpace = (aspRat+1)*tryYSpace;
            
            %Thicken the boundary by lambda plus half rod aspect ratio to ensure no crossing of rods
            %into out-of-bounds area
            outOfBounds = obj.boundaryDesign;
            
            %Calculate lattice, and remove any rods that would fall out of
            %bounds. Then see if number matches target. Loop and reduce
            %lattice spacing until it does. Won't be exactly the density
            %you want, but should be within 3% or so. Can measure more
            %precisely later with obj.getAreaFraction().
            [tryXLat,tryYLat] = meshgrid(1/2*tryXSpace:tryXSpace:obj.xWidth-(1/2*tryXSpace),1/2*tryYSpace:tryYSpace:obj.yHeight-(1/2*tryYSpace));
            inPoints = ~outOfBounds(round(tryYLat(:,1)'*obj.resUp),round(tryXLat(1,:)*obj.resUp));
            
            while sum(inPoints(:)) < tgtRodNo
                tryYSpace = tryYSpace * 0.99;
                tryXSpace = tryXSpace * 0.99;
                
                [tryXLat,tryYLat] = meshgrid(1/2*tryXSpace:tryXSpace:obj.xWidth-(1/2*tryXSpace),1/2*tryYSpace:tryYSpace:obj.yHeight-(1/2*tryYSpace));
                inPoints = ~outOfBounds(round(tryYLat(:,1)'*obj.resUp),round(tryXLat(1,:)*obj.resUp));
            end
            
            %Now actually initialise the rods
            obj.xCells = tryXLat(inPoints);
            obj.yCells = tryYLat(inPoints) + randn(size(obj.xCells))*0.001;
            obj.thetCells = ((rand(size(obj.xCells))>0.5)*pi);
        end
        
        function volFrac = getVolumeFraction(obj) %Calculates the volume fraction occupied by all cells.
            Cylinders = (pi*(obj.lam/2)^2)*(obj.lCells.*(obj.nCells-1));
            Caps = repmat((obj.lam/2)^3 * 4*pi/3,length(obj.nCells),1);
            
            %Need to account for rods forming the barrier - subtract these from the total area estimation
            if ~isempty(obj.lBarr)
                barrCylinders = (pi*(obj.lam/2)^2)*(obj.lBarr.*(obj.nBarr-1));
                barrCaps = repmat((obj.lam/2)^3 * 4*pi/3, length(obj.nBarr),1);
            else
                barrCaps = 0;
                barrCylinders = 0;
            end
            
            totVol = sum(obj.boundaryDesign(:))/(obj.resUp^2);
            volFrac = sum(Caps + Cylinders)/totVol;
        end
        
        function areaFrac = getAreaFraction(obj) %Calculates the area fraction occupied by cells in the z=0 plane. Approximates tilting cells as ellipses. Only works if dz is set to zero.
            %Split barrier and cell populations into those that are above the critical tilt and those below it.
            tiltCritCell = atan(1./(obj.aCells - 1));
            critCells = obj.phiCells > tiltCritCell;
            
            %Calculate subcritical objects as sperocylinders, supercritical objects as ellipses
            arCritCells = (pi * obj.lam^2)./(4*sin(obj.phiCells(critCells)));
            arSubCritCells = (obj.lam^2 * (obj.aCells(~critCells) - 1)) + (pi * (obj.lam/2)^2);
            
            totArea = (obj.xWidth * obj.yHeight) - (sum(obj.boundaryDesign(:))/(obj.resUp^2));
            areaFrac = sum([arCritCells; arSubCritCells])/totArea;
        end
        
        function obj = calcDistMat(obj,contact)
            obj = obj.calcDistThresh(contact);

            if contact
                obj.cellDistsContact = calcGriddedDistMat(obj,strcmp(obj.boundConds,'periodic'));
            else
                obj.cellDistsBiophys = calcGriddedDistMat(obj,strcmp(obj.boundConds,'periodic'));
            end
        end
        
        function obj = calcDistThresh(obj,contact)
            %Calculates the distance threshold below which interactions will be ignored. Based on the longest cell.
            foo = ((obj.nCells - 1) .* obj.lCells) + contact*obj.contactRange;
            maxLen = max(foo); %Length of the longest cell

            obj.distThresh = maxLen + obj.lam + ~contact*log(obj.U0); %Use of log(U0) here is somewhat justified by the exponential drop off in repelling strength. But not terribly. Use cautiously.
        end
        
        function obj = growAndDivide(obj,growthRate,dt,divThresh,postMovement)
            %Grows all cells by a stochastic amount and adds any offspring to the list of exisiting cells.
            obj = obj.growCells(growthRate,dt);
            
            divCells = obj.aCells > divThresh;
            obj = obj.divideCells(divCells);    
            
            %Jiggle colours a little
            obj = obj.jigColours(divCells);
            
            switch postMovement
                case 'same'
                    %Don't actually need to do anything...
                    
                case 'reverse'
                    %Reverse the direction of travel of the original cell (actually the new cell - the one now located at the old end of the cell)
                    obj.thetCells(divCells,1) = mod(obj.thetCells(divCells) + 2*pi,2*pi)-pi;
                    obj.phiCells(divCells,1) = -obj.phiCells(divCells);
                    obj.uCells = [cos(obj.thetCells).*cos(obj.phiCells),sin(obj.thetCells).*cos(obj.phiCells),sin(obj.phiCells)];
                otherwise
                    error('post division movement type not recognized')
            end
        end
        
        function obj = jigColours(obj,divCells)
            %Jiggles the current hue of the object in hue/saturation space and converts back to RGB.
            for i = (size(divCells,1)+1):size(obj.cCells,1)
                outmap = rgb2hsv(obj.cCells(i,:));
                
                outmap(1) = mod(outmap(1) + (randn(1)*obj.colJigRate),1);
                obj.cCells(i,:) = hsv2rgb(outmap);
            end
        end
        
        function obj = setColours(obj,colourCells)
            %Sets the colours of the cells according to position.
            %Do any cell colouration steps here
            switch colourCells
                case 'Tilt'
                    obj.cCells = [obj.phiCells/pi,obj.phiCells/pi,1-(obj.phiCells/pi)]; %Tilt
                case 'Orientation'
                    map = colormap('hsv');
                    newColors = zeros(size(obj.cCells));
                    for k = 1:size(obj.thetCells,1)
                        bin = ceil((mod(obj.thetCells(k)+pi,pi)/pi) * size(map,1));
                        newColors(k,:) = map(bin,:);
                    end
                    obj.cCells = newColors;
                case 'Position'
                    for k = 1:size(obj.cCells,1)
                        xFac = obj.xCells(k)/obj.xWidth;
                        yFac = obj.yCells(k)/obj.yHeight;
                        obj.cCells(k,:) = [xFac,1-(xFac+yFac)/2,yFac];
                    end
                case 'Population'
                    popArr = [obj.aCells,obj.fCells,obj.fireCells,obj.rCells]; %Believe these to be the only properties that could distinguish populations
                    popArr = sortrows(popArr);
                    [~,~,ic] = unique(popArr,'rows');
                    cmap = colormap('lines');
                    ic(ic > size(cmap,1)) = size(cmap,1);
                    obj.cCells = cmap(ic,:);
                case 'Hits' %Choice for debugging purposes - won't necessarily always work with all parameter settings. Assumes a firing and non-firing population.
                    maxHit = 5;
                    for k = 1:size(obj.cCells,1)
                        currHits = min(obj.hitCells(k),maxHit);
                        if currHits > 0 %Dead cells are set as magenta shading to gray
                            colFac = currHits/(maxHit+0.5);
                            obj.cCells(k,:) = 1-[0.9,colFac*0.9,0.4+colFac*0.5];
                        elseif obj.fireCells(k) > 0 %Firing cells are set as light blue
                            obj.cCells(k,:) = 1-([32,128,196]/255);
                        elseif obj.hitCells(k) == 0 %Unhit cells are set as yellow
                            obj.cCells(k,:) = 1-([255, 190, 11]/255);
                        end
                    end
            end
        end
        
        function obj = divideCells(obj,divInds)
            if sum(divInds) > 0
                xFac = cos(obj.thetCells(divInds)) .* cos(obj.phiCells(divInds)) .* obj.lCells(divInds) .* (obj.nCells(divInds) + 1) ./ 4; %To get the separation, should really subtract 1 from nCells. But this leads to unpleasantly close cells.
                yFac = sin(obj.thetCells(divInds)) .* cos(obj.phiCells(divInds)) .* obj.lCells(divInds) .* (obj.nCells(divInds) + 1) ./ 4;
%                 zFac = cos(obj.phiCells(divInds)) .* obj.lCells(divInds) .* (obj.nCells(divInds) + 1) ./ 4;
                oldX = obj.xCells(divInds) - xFac;
                oldY = obj.yCells(divInds) - yFac;
%                 obj.zCells(divInds) = obj.zCells(divInds) - zFac;
                obj.aCells(divInds) = obj.aCells(divInds)/2;
                obj.hitCells(divInds) = obj.hitCells(divInds)/2;
                
                %Add in new cells
                newX = oldX + 2*xFac;
                newY = oldY + 2*yFac;
                newZ = obj.zCells(divInds);
                newA = obj.aCells(divInds);
                newL = obj.lCells(divInds);
                newThet = obj.thetCells(divInds) + (randn(sum(divInds),1)*0.2);
                newPhi = -obj.phiCells(divInds);
                newN = obj.nCells(divInds);
                newU = obj.uCells(divInds,:);
                newC = obj.cCells(divInds,:);
                newF = obj.fCells(divInds,:);
                newFire = obj.fireCells(divInds,:);
                newHit = obj.hitCells(divInds,:);
                
                %Need to deal with periodic boundary conditions
                if strcmp(obj.boundConds,'periodic')
                    oldX(oldX < 0) = oldX(oldX < 0) + obj.xWidth;
                    oldY(oldY < 0) = oldY(oldY < 0) + obj.yHeight;
                    
                    oldX(oldX > obj.xWidth) = oldX(oldX > obj.xWidth) - obj.xWidth;
                    oldY(oldY > obj.yHeight) = oldY(oldY > obj.yHeight) - obj.yHeight;
                    
                    newX(newX < 0) = newX(newX < 0) + obj.xWidth;
                    newY(newY < 0) = newY(newY < 0) + obj.yHeight;
                    
                    newX(newX > obj.xWidth) = newX(newX > obj.xWidth) - obj.xWidth;
                    newY(newY > obj.yHeight) = newY(newY > obj.yHeight) - obj.yHeight;
                end
                
                obj.xCells(divInds) = oldX;
                obj.yCells(divInds) = oldY;
                
                obj.xCells = [obj.xCells;newX];
                obj.yCells = [obj.yCells;newY];
                obj.zCells = [obj.zCells;newZ];
                obj.aCells = [obj.aCells;newA];
                obj.lCells = [obj.lCells;newL];
                obj.thetCells = [obj.thetCells;newThet];
                obj.phiCells = [obj.phiCells;newPhi];
                obj.nCells = [obj.nCells;newN];
                obj.uCells = [obj.uCells;newU];
                obj.cCells = [obj.cCells;newC];
                obj.fCells = [obj.fCells;newF];
                obj.hitCells = [obj.hitCells;newHit];
                obj.fireCells = [obj.fireCels;newFire];
                
                %Recalculate segment properties.
                [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
                %             obj = obj.calculateSegmentPositions();
            end
        end
        
        function obj = growCells(obj,growthRate,dt)
            %Adds length to the cell at the rate specified by growthRate. Growth is proportional to cell length.
            addLength = (growthRate*dt)*abs(randn(size(obj.aCells))).*obj.aCells;
            obj.aCells = obj.aCells + addLength;
            
            [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
            
            %Update segment positions (introducing new segments if needed).
%             obj = obj.calculateSegmentPositions();
        end
        
        function contacts = calculateContacts(obj)
            %Calculates a list of cells that are touching each other, as
            %defined in the CDI hit detection code.
%             obj = obj.calcDistMat(true);
% 
%             contactMat = zeros(size(obj.xCells,1));
%             for i = 1:size(obj.xCells,1)
%                 for j = (i+1):size(obj.xCells,1)
%                     if ~isnan(obj.cellDistsContact(i,j))
%                         contactMat(i,j) = testEllipseCollision(obj,i,j);
%                         contactMat(j,i) = contactMat(i,j);
%                     end
%                 end
%             end
% 
%             contacts = cell(size(obj.xCells,1),1);
%             for i = 1:size(obj.xCells,1)
%                 contacts{i} = find(contactMat(:,i));
%             end
%             %Calculates a list of cells that are touching each other, as
%             %defined in the CDI hit detection code.
            pxSize = 0.25; %Granularity of the pixel-based approach for calculating elliptical neighbourhoods. Set smaller to improve resolution of approach.
            contacts = cell(size(obj.xCells,1),1);
            
            %Begin by creating an image of the indices of each cell
            indexImg = paintOverlapField(obj,pxSize);
            
            %Now step through each rod, and find which other rods it is
            %within CDI range of. Apply hits to these probabilistically,
            %based on the focal rod's firing rate.
            for i = 1:size(obj.xCells,1)
                expandImg = makeExpandedRodProfile(obj,i,pxSize);
                
                contactInds = unique(indexImg(expandImg));
                contactInds(contactInds == 0) = [];
                contactInds(contactInds == i) = [];
                
                contacts{i} = contactInds;
            end
        end
        
        function [obj,hitNoNew,hitNoTot] = calculateHits(obj,dt)
            %Calculates hits from a contact-dependant killing system (e.g.
            %CDI, T6SS). Note - assumes sampling rate (dt) is much faster
            %than the hit rate, so a maximum of one hit per timepoint (per
            %cell-cell interaction) is applied.
            pxSize = 0.25; %Granularity of the pixel-based approach for calculating elliptical neighbourhoods. Set smaller to improve resolution of approach.

            %Begin by creating an image of the indices of each cell
            indexImg = paintOverlapField(obj,pxSize);

            %Now step through each rod, and find which other rods it is
            %within CDI range of. Apply hits to these probabilistically,
            %based on the focal rod's firing rate.
            hitNoNew = 0;
            hitNoTot = 0;
            for i = 1:size(obj.xCells,1)
                if obj.fireCells(i) > 0 %If this is an attacker
                    expandImg = makeExpandedRodProfile(obj,i,pxSize);

                    contactInds = unique(indexImg(expandImg));
                    contactInds(contactInds == 0) = [];
                    contactInds(contactInds == i) = [];

%             obj = obj.calcDistMat(true);
% 
%             for i = 1:size(obj.xCells,1)
%                 if obj.fireCells(i) > 0 %If this is an attacker
%                     
%                     contactInds = [];
%                     for j = 1:size(obj.xCells)
%                         if i ~= j && ~isnan(obj.cellDistsContact)
%                             if testEllipseCollision(obj,i,j)
%                                 contactInds = [contactInds;j];
%                             end
%                         end
%                     end

                    switch obj.hitRateType
                        case 'distributed'
                            hitProb = obj.fireCells(i)*dt/numel(contactInds);
                        case 'constant'
                            hitProb = obj.fireCells(i)*dt;
                    end

                    hitEvts = rand(size(contactInds)) < hitProb;
                    hitInds = contactInds(hitEvts);

                    %This bit of code prevents you from accumulating hits
                    %from cells of the same population - effectively
                    %simulating cognate immunity gene expression
                    hitPops = obj.popCells(hitInds);
                    thisPop = obj.popCells(i);
                    hitInds = hitInds(hitPops ~= thisPop);

                    hitNoNew = hitNoNew + sum(obj.hitCells(hitInds) == 0); %Total number of cells that have never been hit that are about to be hit for the first time
                    hitNoTot = hitNoTot + size(hitInds,1); %Total number of hits applied by this attacker cell to sensitive cells

                    obj.hitCells(hitInds) = obj.hitCells(hitInds) + 1;
                end
            end
        end
        
        function obj = killCells(obj)
            %Kills cells on the basis of the number of hits from the toxin
            %system they have recieved. Any cells that exceed the specified
            %threshold are removed from the simulation (killType = lyse)
            %or become inactive 'husks' (killType = husk)
            killInds = obj.hitCells >= obj.killThresh;

            switch obj.killType
                case 'lyse'
                    %Remove killed cells from all fields
                    obj.xCells(killInds) = [];
                    obj.yCells(killInds) = [];
                    obj.zCells(killInds) = [];
                    obj.thetCells(killInds) = [];
                    obj.phiCells(killInds) = [];
                    obj.aCells(killInds) = [];
                    obj.nCells(killInds) = [];
                    obj.uCells(killInds,:) = [];
                    obj.lCells(killInds) = [];
                    obj.fCells(killInds) = [];
                    obj.rCells(killInds) = [];
                    obj.cCells(killInds,:) = [];
                    obj.hitCells(killInds) = [];
                    obj.fireCells(killInds) = [];
                    obj.popCells(killInds) = [];
                    
                    obj = obj.calcDistMat(false);
                case 'husk'
                    %Inactivate all active behaviours of killed cells
                    obj.fCells(killInds) = 0;
                    obj.rCells(killInds) = 0;
                    obj.fireCells(killInds) = 0;
                    obj.popCells(killInds) = 'd';
            end
        end
        
        function [drdt,dthetadt,dphidt] = calcVelocities(obj,stepType)
            %Calculates the rate of translation and rotation for all cells             
            [fT,fR,fPar] = calcFrictionTensors(obj.aCells,obj.uCells,obj.f0,obj.frictionAnisotropy);
            
            %I will provide three methods for calculating the potential
            %between rods (the slowest part of the model). Option 1, the
            %most speedy, is to use the graphics card to split the
            %calculation between many workers. Option 2, less speedy, is
            %to use a pre-compiled .mex file. Option 3, less speedy still,
            %is to use Matlab's inbuilt functions. Try each option in turn,
            %opting for the next if the necessary hardware/code is not
            %available.
            if obj.gpuAvailable
                [gradXYZ,gradTheta,gradPhi] = calcPotentialGradsGPU(obj);
            elseif obj.compiled
                [gradXYZ,gradTheta,gradPhi] = calcPotentialGradsCompiled(obj);
            else
                [gradXYZ,gradTheta,gradPhi] = calcPotentialGradsBase(obj);
            end
            
            drdt = zeros(size(obj.xCells,1),3);
            dthetadt = zeros(size(obj.xCells));
            dphidt = zeros(size(obj.xCells));
                
            for i = 1:size(obj.xCells,1)
                v0 = obj.fCells(i)/(obj.f0*fPar(i)); %The self-propulsion velocity of a non-interacting SPR
                uAlph = obj.uCells(i,:)'; %The unit vector representing the current orientation of the rod
                
                drdt(i,:) = (v0*uAlph - (fT(:,:,i)\gradXYZ(i,:)'))';
                drdt(i,3) = 0; %Set the z movement component to be 0 (so cells don't move out of the monolayer)
                dthetadt(i) = -gradTheta(i)/fR(i);
                if ~isinf(obj.zElasticity) && ~strcmp(stepType,'burnIn') && ~strcmp(stepType,'Initial')
                    dphidt(i) = -(gradPhi(i) + (obj.zElasticity/2) * (obj.aCells(i)^2) * cos(obj.phiCells(i)) * sin(obj.phiCells(i)))/fR(i); %We use a Kelvin-Voigt model, with gradPhi being taken as proportional to the stress term and strain proportional to the displacement of the cell from the neutral (phi = 0) position
                else
                    dphidt(i) = 0;
                end
            end
        end
        
        function obj = stepModel(obj,timeStep,growthRate,divThresh,postMovement,colourCells,stepType)
            %Increases the time by one step, updating cell positions based on their current velocities.
            
            %Calculate distance matrix and threshold if needed (e.g. for first time point)
            if isempty(obj.cellDistsBiophys)
                obj = obj.calcDistMat(false);
            end
               
            %Need to re-calculate distance threshold at each step. May change as cells grow.
            
%             %Apply the RK4 method to simulate movement dynamics
%             [drdtk1,dthetadtk1,dphidtk1] = obj.calcVelocities(stepType);
%             k2 = obj.moveCells(drdtk1,dthetadtk1,dphidtk1,timeStep/2);
%             [drdtk2,dthetadtk2,dphidtk2] = k2.calcVelocities(stepType);
%             k3 = obj.moveCells(drdtk2,dthetadtk2,dphidtk2,timeStep/2);
%             [drdtk3,dthetadtk3,dphidtk3] = k3.calcVelocities(stepType);
%             k4 = obj.moveCells(drdtk3,dthetadtk3,dphidtk3,timeStep);
%             [drdtk4,dthetadtk4,dphidtk4] = k4.calcVelocities(stepType);
%             
%             drdt = 1/6*(drdtk1 + 2*drdtk2 + 2*drdtk3 + drdtk4);
%             dthetadt = 1/6*(dthetadtk1 + 2*dthetadtk2 + 2*dthetadtk3 + dthetadtk4);
%             dphidt = 1/6*(dphidtk1 + 2*dphidtk2 + 2*dphidtk3 + dphidtk4);
            
            %Apply midpoint method to simulate movement dynamics
            [drdtk1,dthetadtk1,dphidtk1] = obj.calcVelocities(stepType);
            k1 = obj.moveCells(drdtk1,dthetadtk1,dphidtk1,timeStep/2);
            [drdt,dthetadt,dphidt] = k1.calcVelocities(stepType);
            
%             %Apply Euler method to simulate movement dynamics
%             [drdt,dthetadt,dphidt] = obj.calcVelocities(stepType);

            %Update position and angle of cells based on dynamic parameters and timestep size
            obj = obj.moveCells(drdt,dthetadt,dphidt,timeStep);
            
            %Update the lengths of the cells and add any extra daughter cells that arise.
            obj = obj.growAndDivide(growthRate,timeStep,divThresh,postMovement);
            
            %Invert the movement direction of any cells regarded to have reversed.
            obj = obj.randomReverse(timeStep);
            
            obj = obj.setColours(colourCells);
            
            %Calculate distance matrix. Allows elimination of small-intensity interactions at later time points.
            obj = obj.calcDistMat(false);
        end
        
        function obj = moveCells(obj,drdt,dthetadt,dphidt,timeStep)
            %Updates the position of the cell based on current velocity
            obj.xCells = obj.xCells + timeStep*drdt(:,1);
            obj.yCells = obj.yCells + timeStep*drdt(:,2);
            obj.zCells = obj.zCells + timeStep*drdt(:,3);
            
            if strcmp(obj.boundConds,'periodic')
                %Apply periodic boundary conditions
                obj.xCells = mod(obj.xCells,obj.xWidth);
                obj.yCells = mod(obj.yCells,obj.yHeight);
                obj.zCells = mod(obj.zCells,obj.zDepth);
            end
            
            %Updates the angles of the cell
            obj.thetCells = obj.thetCells + dthetadt*timeStep;
            obj.thetCells = mod(obj.thetCells + pi,2*pi) - pi;
            obj.phiCells = obj.phiCells + dphidt*timeStep;
            obj.phiCells = mod(obj.phiCells + pi/2, pi) - pi/2;
            
            obj.uCells = [cos(obj.thetCells).*cos(obj.phiCells),sin(obj.thetCells).*cos(obj.phiCells),sin(obj.phiCells)];
        end
        
        function obj = randomReverse(obj,timeStep)
            %Reverse rate is the probability of a given cell reversing in a single time step. Results in a Binomial distribution of reversal events for a given cell (?).
            reversingCells = obj.rCells > rand(size(obj.rCells))/timeStep;
            obj.thetCells(reversingCells) = rem(obj.thetCells(reversingCells) + 2*pi,2*pi) - pi;
        end
        
        function outImg = drawField(obj)
            %Draws the current state of the model - location and angles of all rods in model. Coloured rods are motile cells.
            %Note that this is only a 2D projection for debugging. For proper imaging of the 3D system, use the paraview scripts.
            Upsample = 20; %Extent to which the 'design' image should be interpolated to create smoother graphics.
            Downsample = 5; %Extent to which the final image should be scaled down to save space.
            fullSF = Upsample * obj.resUp; %Extent to which the ellipse images should be blown up.

            majorBoost = 0.2; %Extra factor by which to make cells longer for rendering purposes
            minorBoost = 0.1; %Extra factor by which to make cells wider for rendering purposes
            
%             hand = figure(1);
%             set(hand, 'Units', 'pixels', 'Position', posVec);
%             axis([0,obj.xWidth,0,obj.yHeight])
%             hand.CurrentAxes.Position = [0,0,1,1];
%             hand.CurrentAxes.XTick = [];
%             hand.CurrentAxes.YTick = [];
%             set(gca,'YDir','reverse')
%             colorbar off
%             cla
%             hold on
            
            %Draw barrier (as an image) and break into separate rgb
            %channels
            background = imresize(double(~obj.boundaryDesign),Upsample); %We'll do the smoothing later
            imgr = imgaussfilt(double(background),fullSF*obj.lam);
            imgg = imgaussfilt(double(background),fullSF*obj.lam);
            imgb = imgaussfilt(double(background),fullSF*obj.lam);
            
            %Draw rods
            
            %Do separate paintjobs for each of the three colour
            %channels
            imgr = paintEllipse(imgr,obj.xCells,obj.yCells,(obj.aCells/2)+majorBoost,(obj.lam*ones(size(obj.aCells))/2)+minorBoost,-rad2deg(obj.thetCells),obj.cCells(:,1),1/fullSF,obj.boundConds,obj.xWidth,obj.yHeight);
            imgg = paintEllipse(imgg,obj.xCells,obj.yCells,(obj.aCells/2)+majorBoost,(obj.lam*ones(size(obj.aCells))/2)+minorBoost,-rad2deg(obj.thetCells),obj.cCells(:,2),1/fullSF,obj.boundConds,obj.xWidth,obj.yHeight);
            imgb = paintEllipse(imgb,obj.xCells,obj.yCells,(obj.aCells/2)+majorBoost,(obj.lam*ones(size(obj.aCells))/2)+minorBoost,-rad2deg(obj.thetCells),obj.cCells(:,3),1/fullSF,obj.boundConds,obj.xWidth,obj.yHeight);
            
            %Smooth to make nicer looking
            imgr = imgaussfilt(imgr,fullSF*obj.lam/4);
            imgg = imgaussfilt(imgg,fullSF*obj.lam/4);
            imgb = imgaussfilt(imgb,fullSF*obj.lam/4);
            
            outImg = cat(3,imgr,imgg,imgb);
%             imshow(outImg,'Parent',hand.CurrentAxes);
%             axis(hand.CurrentAxes,'tight')
            
            outImg = imresize(outImg,1/Downsample);
        end
        
        function axHand = plotField(obj,posVec,drawContacts)
            %Draws the current state of the model - location and angles of all cells in model. Black spots are front of cells.
            %Note - this is the old (inelegant) method which is more
            %compatible with Matlab's other plotting functions. For nice,
            %quick visuals, use the drawField function
            figHand = figure(1);
            set(figHand, 'Units', 'pixels', 'Position',posVec);
            axis([0,obj.xWidth,0,obj.yHeight])
            axHand = figHand.CurrentAxes;
            axHand.Position = [0,0,1,1];
            axHand.XTick = [];
            axHand.YTick = [];
            set(gca,'YDir','reverse')
            colorbar off
            cla
            hold on
            plot(obj.xWidth-(obj.xWidth/posVec(3)),obj.yHeight-(obj.yHeight/posVec(4)),'r.')
            
            for i = 1:length(obj.xCells)
                h = ellipse((obj.lCells(i)*(obj.nCells(i) - 1) + obj.lam)*0.5,obj.lam*0.5,obj.thetCells(i),obj.xCells(i),obj.yCells(i),[0,0,0]);
                x=get(h,'Xdata');
                y=get(h,'Ydata');
                delete(h)
                patch(x,y,1-obj.cCells(i,:),'EdgeColor','none');
                
                if strcmp(obj.boundConds,'periodic')
                    %Also draw an ellipse at each of the periodic locations. Won't be rendered unless within field of view.
                    if obj.xCells(i) - (obj.aCells(i)*obj.lam) < 0
                        h = ellipse((obj.lCells(i)*(obj.nCells(i) - 1) + obj.lam)*0.5,obj.lam*0.5,obj.thetCells(i),obj.xCells(i) + obj.xWidth,obj.yCells(i),[0,0,0]);
                        x=get(h,'Xdata');
                        y=get(h,'Ydata');
                        delete(h)
                        patch(x,y,1-obj.cCells(i,:),'EdgeColor','none');
                    elseif obj.xCells(i) + (obj.aCells(i)*obj.lam) > obj.xWidth
                        h = ellipse((obj.lCells(i)*(obj.nCells(i) - 1) + obj.lam)*0.5,obj.lam*0.5,obj.thetCells(i),obj.xCells(i) - obj.xWidth,obj.yCells(i),[0,0,0]);
                        x=get(h,'Xdata');
                        y=get(h,'Ydata');
                        delete(h)
                        patch(x,y,1-obj.cCells(i,:),'EdgeColor','none');
                    end
                    
                    if obj.yCells(i) - (obj.aCells(i)*obj.lam) < 0
                        h = ellipse((obj.lCells(i)*(obj.nCells(i) - 1) + obj.lam)*0.5,obj.lam*0.5,obj.thetCells(i),obj.xCells(i),obj.yCells(i) + obj.yHeight,[0,0,0]);
                        x=get(h,'Xdata');
                        y=get(h,'Ydata');
                        delete(h)
                        patch(x,y,1-obj.cCells(i,:),'EdgeColor','none');
                    elseif obj.yCells(i) + (obj.aCells(i)*obj.lam) > obj.yHeight
                        h = ellipse((obj.lCells(i)*(obj.nCells(i) - 1) + obj.lam)*0.5,obj.lam*0.5,obj.thetCells(i),obj.xCells(i),obj.yCells(i) - obj.yHeight,[0,0,0]);
                        x=get(h,'Xdata');
                        y=get(h,'Ydata');
                        delete(h)
                        patch(x,y,1-obj.cCells(i,:),'EdgeColor','none');
                    end
                end
            end
            
            if drawContacts
                contacts = obj.calculateContacts();
                for i = 1:size(contacts,1)
                    for j = 1:size(contacts{i},1)
                        x1 = obj.xCells(i);
                        y1 = obj.yCells(i);
                        y2 = obj.yCells(contacts{i}(j));
                        x2 = obj.xCells(contacts{i}(j));

                        if abs(x1-x2) > obj.xWidth/2
                            if x1 < x2
                                x1 = x1 + obj.xWidth;
                            else
                                x2 = x2 + obj.xWidth;
                            end
                        end

                        if abs(y1-y2) > obj.yHeight/2
                            if y1 < y2
                                y1 = y1 + obj.yHeight;
                            else
                                y2 = y2 + obj.yHeight;
                            end
                        end

                        plot(x1,y1,'w.','MarkerSize',15)
                        plot(x2,y2,'w.','MarkerSize',15)
                        plot([x1,x2],[y1,y2],'k','LineWidth',1.5)
                    end
                end
            end
        end
        
    end
end