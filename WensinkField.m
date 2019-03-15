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
        
        xBarr %x position for the barrier rods (non-motile, but still cause steric effects)
        yBarr
        zBarr
        thetBarr
        phiBarr
        aBarr %The approach here will be to model the boundary segments as rods of aspect ratio a = 1. 
        nBarr %This is not the most elegant solution, but allows us to recycle the code for rod-rod interactions more easily.
        lBarr
        
        cellDists %Distance between all cells. Indexed in same way as obj.cells.
        distThresh %Distance between cells, beyond which interactions are too weak to be relevant. 
        U0 %Potential Amplitude
        lam %Screening length of Yukawa segements
        boundConds %Boundary conditions (periodic or none)
        boundaryDesign %Image containing an image of the design of the confinement
        resUp %Difference in resolution between boundaryDesign image and actual simulation domain
    end
    methods
        function obj = WensinkField(xWidth,yHeight,zDepth,U0,lam,boundaryConditions)
            if isnumeric(xWidth) && xWidth > 0 && isnumeric(yHeight) && yHeight > 0
                obj.xWidth = xWidth;
                obj.yHeight= yHeight;
                obj.zDepth = zDepth;
            else
                error('Input arguments to WensinkField are not valid');
            end
            
            obj.U0 = U0;
            obj.lam = lam;
            obj.boundConds = boundaryConditions;
        end
        
        function obj = populateField(obj,barrierSettingsType,barrierSettings,cellSettingsType,cellSettings,areaFrac)
            %Populates the field with a number of randomly positioned cells of the given force generation and a number of static 'barrier' rods.
            
            switch barrierSettingsType
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
            
            switch cellSettingsType
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
                case 'LatticedXYCells' %Start with an initially ordered lattice of cells (with random up/down orientations)
                    obj = obj.overlayLattice(areaFrac,cellSettings.a);
                    obj.zCells = ones(size(obj.xCells)) * obj.zDepth/2;
                    obj.phiCells = abs(randn(size(obj.xCells)))*0.001;
                    obj.aCells = cellSettings.a * ones(size(obj.xCells));
                    obj.fCells = cellSettings.f * ones(size(obj.xCells));
                    obj.rCells = cellSettings.r * ones(size(obj.xCells));
                    obj.cCells = rand(size(obj.xCells,1),3);
              end
            
            [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
            obj.uCells = [cos(obj.thetCells).*cos(obj.phiCells),sin(obj.thetCells).*cos(obj.phiCells),sin(obj.phiCells)];
        end
        
        function obj = overlayLattice(obj,areaFrac,aspRat)
            %Overlays a lattice of rods on top of the region defined by the
            %design image, and adjusts density until it matches the desired
            %area fraction. Assumes all rods are identical in size.
            totArea = sum(~obj.boundaryDesign(:))/(obj.resUp^2);
            singleArea = (obj.lam^2 * (aspRat - 1)) + (pi * (obj.lam/2)^2);
            
            tgtRodNo = round(totArea*areaFrac/singleArea);
            
            %Approach will be to do a lattice that is slightly too loose,
            %and then gradually tighten it until the target no. of rods fit within
            %the bounded domain.
            tryYSpace = sqrt(singleArea/(areaFrac*(aspRat+1))); %Spacing of lattice in y-direction
            tryXSpace = (aspRat+1)*tryYSpace;
            
            %Thicken the boundary by lambda plus half rod aspect ratio to ensure no crossing of rods
            %into out-of-bounds area
            outOfBounds = bwmorph(obj.boundaryDesign,'thicken',round(obj.resUp*(obj.lam+(aspRat/2))));
            
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
            
            totArea = sum(~obj.boundaryDesign(:))/(obj.resUp^2);
            areaFrac = sum([arCritCells; arSubCritCells])/totArea;
        end
        
        function obj = calcDistMat(obj)
            if strcmp(obj.boundConds,'periodic')
                obj.cellDists = calcGriddedDistMat(obj,true);
            else
                obj.cellDists = calcGriddedDistMat(obj,false);
            end
        end
        
        function obj = calcDistThresh(obj)
            %Calculates the distance threshold below which interactions will be ignored. Based on the longest cell.
            foo = ((obj.nCells - 1) .* obj.lCells);
            maxLen = max(foo); %Length of the longest cell
            
            obj.distThresh = maxLen + obj.lam + log(obj.U0); %Use of log(U0) here is somewhat justified by the exponential drop off in repelling strength. But not terribly. Use cautiously.
        end
        
        function obj = growAndDivide(obj,growthRate,dt,divThresh,postMovement,colJigRate)
            %Grows all cells by a stochastic amount and adds any offspring to the list of exisiting cells.
            obj = obj.growCells(growthRate,dt);
            
            divCells = obj.aCells > divThresh;
            obj = obj.divideCells(divCells);    
            
            %Jiggle colours a little
            obj = obj.jigColours(colJigRate,divCells);
            
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
        
        function obj = jigColours(obj,jigRate,divCells)
            %Jiggles the current hue of the object in hue/saturation space and converts back to RGB.
            for i = (size(divCells,1)+1):size(obj.cCells,1)
                outmap = rgb2hsv(obj.cCells(i,:));
                
                outmap(1) = mod(outmap(1) + (randn(1)*jigRate),1);
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
                        obj.cCells(k,:) = [xFac,1,yFac];
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
        
        function [drdt,dthetadt,dphidt] = calcVelocities(obj,f0,zElasticity)
            %Calculates the rate of translation and rotation for all cells
            includeMat = obj.cellDists < obj.distThresh;
            includeMat(logical(diag(ones(length(includeMat),1)))) = 0;
            includeMat = includeMat(1:length(obj.nCells),:); %Alpha rods should include all cells, beta rods all cells and all barrier rods
             
            [fT,fR] = calcFrictionTensors(obj.aCells,obj.uCells,f0);
            [fPar,~,~] = calcGeomFactors(obj.aCells);
            
            drdt = zeros(size(obj.xCells,1),3);
            dthetadt = zeros(size(obj.xCells));
            dphidt = zeros(size(obj.xCells));
            
            boundX = obj.xWidth/2;
            boundY = obj.yHeight/2;
            boundZ = obj.zDepth/2;
            Height = obj.yHeight;
            Width = obj.xWidth;
            Depth = obj.zDepth;
            lam = obj.lam;
            U0 = obj.U0;
            
            for i = 1:length(obj.nCells) %The index of the cell alpha.
                %Get indices of cells that this cell (alpha) interacts with
                betInds = find(includeMat(i,:));
                xs = [obj.xCells;obj.xBarr]; xBets = xs(betInds);
                ys = [obj.yCells;obj.yBarr]; yBets = ys(betInds);
                zs = [obj.zCells;obj.zBarr]; zBets = zs(betInds);
                ns = [obj.nCells;obj.nBarr]; nBets = ns(betInds);
                ls = [obj.lCells;obj.lBarr]; lBets = ls(betInds);
                thets = [obj.thetCells;obj.thetBarr]; thetBets = thets(betInds);
                phis = [obj.phiCells;obj.phiBarr]; phiBets = phis(betInds);
                
                %Get dynamics

                if length(obj.aCells) >= 1
                    xAlph = obj.xCells(i);
                    yAlph = obj.yCells(i);
                    zAlph = obj.zCells(i);
                    lAlph = obj.lCells(i);
                    nAlph = obj.nCells(i);
                    thetAlph = obj.thetCells(i);
                    uAlph = obj.uCells(i,:)';
                    phiAlph = obj.phiCells(i);
                end
                
                if strcmp(obj.boundConds,'none')
                    [dUdx,dUdy,dUdz,dUdthet,dUdphi] = mexCalcEnergyGradients(xBets,yBets,zBets,lBets,nBets,thetBets,phiBets,xAlph,yAlph,zAlph,lAlph,nAlph,thetAlph,phiAlph,U0,lam,boundX,boundY,Width,Height);
                elseif strcmp(obj.boundConds,'periodic')
                    [dUdx,dUdy,dUdz,dUdthet,dUdphi] = mexCalcEnergyGradientsPeriodic(xBets,yBets,zBets,lBets,nBets,thetBets,phiBets,xAlph,yAlph,zAlph,lAlph,nAlph,thetAlph,phiAlph,U0,lam,boundX,boundY,Width,Height);
                end
                
                if sum(isnan(dUdx)) > 0
                    disp('break')
                end
                
                gradXY = -[sum(dUdx)/2,sum(dUdy)/2,sum(dUdz)/2];
                gradTheta = -sum(dUdthet)/2;
                gradPhi = -sum(dUdphi)/2;
                
                v0 = obj.fCells(i)/(f0*fPar(i)); %The self-propulsion velocity of a non-interacting SPR
                
                drdt(i,:) = (v0*uAlph - (fT(:,:,i)\gradXY'))';
                drdt(i,3) = 0; %Set the z movement component to be 0 (so cells don't move out of the monolayer)
                dthetadt(i) = -gradTheta/fR(i);
                if ~isinf(zElasticity) 
                    dphidt(i) = -(gradPhi + (zElasticity/2) * (obj.aCells(i)^2) * cos(phiAlph) * sin(phiAlph))/fR(i); %We use a Kelvin-Voigt model, with gradPhi being taken as proportional to the stress term and strain proportional to the displacement of the cell from the neutral (phi = 0) position
                else
                    dphidt(i) = 0;
                end
            end
        end
        
        function obj = stepModel(obj,timeStep,f0,zElasticity,growthRate,divThresh,postMovement,colJigRate,colourCells)
            %Increases the time by one step, updating cell positions based on their current velocities.
            
            %Calculate distance matrix and threshold if needed (e.g. for first time point)
            if isempty(obj.distThresh)
                obj = obj.calcDistThresh();
            end
            if isempty(obj.cellDists)
                obj = obj.calcDistMat();
            end
               
            %Need to re-calculate distance threshold at each step. May change as cells grow.
            
%             %Apply the RK4 method to simulate movement dynamics
%             [drdtk1,dthetadtk1,dphidtk1] = obj.calcVelocities(f0,zElasticity);
%             k2 = obj.moveCells(drdtk1,dthetadtk1,dphidtk1,timeStep/2);
%             [drdtk2,dthetadtk2,dphidtk2] = k2.calcVelocities(f0,zElasticity);
%             k3 = obj.moveCells(drdtk2,dthetadtk2,dphidtk2,timeStep/2);
%             [drdtk3,dthetadtk3,dphidtk3] = k3.calcVelocities(f0,zElasticity);
%             k4 = obj.moveCells(drdtk3,dthetadtk3,dphidtk3,timeStep);
%             [drdtk4,dthetadtk4,dphidtk4] = k4.calcVelocities(f0,zElasticity);
%             
%             drdt = 1/6*(drdtk1 + 2*drdtk2 + 2*drdtk3 + drdtk4);
%             dthetadt = 1/6*(dthetadtk1 + 2*dthetadtk2 + 2*dthetadtk3 + dthetadtk4);
%             dphidt = 1/6*(dphidtk1 + 2*dphidtk2 + 2*dphidtk3 + dphidtk4);
            
            %Apply midpoint method to simulate movement dynamics
            [drdtk1,dthetadtk1,dphidtk1] = obj.calcVelocities(f0,zElasticity);
            k1 = obj.moveCells(drdtk1,dthetadtk1,dphidtk1,timeStep/2);
            [drdt,dthetadt,dphidt] = k1.calcVelocities(f0,zElasticity);
            
%             %Apply Euler method to simulate movement dynamics
%             [drdt,dthetadt,dphidt] = obj.calcVelocities(f0,zElasticity);

            %Update position and angle of cells based on dynamic parameters and timestep size
            obj = obj.moveCells(drdt,dthetadt,dphidt,timeStep);
            
            tmp = obj.growAndDivide(growthRate,timeStep,divThresh,postMovement,colJigRate);
            if sum(tmp.xCells < 0) > 0 || sum(tmp.yCells < 0) > 0 || sum(tmp.xCells > tmp.xWidth) > 0 || sum(tmp.yCells > tmp.yHeight) > 0
                disp('break')
            end
            
            %Update the lengths of the cells and add any extra daughter cells that arise.
            obj = obj.growAndDivide(growthRate,timeStep,divThresh,postMovement,colJigRate);
            
            if sum(obj.xCells < 0) > 0 || sum(obj.yCells < 0) > 0 || sum(obj.xCells > obj.xWidth) > 0 || sum(obj.yCells > obj.yHeight) > 0
                disp('break')
            end
            
            %Invert the movement direction of any cells regarded to have reversed.
            obj = obj.randomReverse();
            
            obj = obj.setColours(colourCells);
            
            %Calculate distance matrix. Allows elimination of small-intensity interactions at later time points.
            obj = obj.calcDistThresh();
            obj = obj.calcDistMat();
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
        
        function obj = randomReverse(obj)
            %Reverse rate is the probability of a given cell reversing in a single time step. Results in a Binomial distribution of reversal events for a given cell (?).
            reversingCells = obj.rCells > rand(size(obj.rCells));
            obj.thetCells(reversingCells) = rem(obj.thetCells(reversingCells) + 2*pi,2*pi) - pi;
        end
        
        function outImg = drawField(obj,posVec)
            %Draws the current state of the model - location and angles of all rods in model. Coloured rods are motile cells.
            %Note that this is only a 2D projection for debugging. For proper imaging of the 3D system, use the paraview scripts.
            Upsample = 5; %Extent to which the 'design' image should be interpolated to create smoother graphics.
            Downsample = 2.5; %Extent to which the final image should be scaled down to save space.
            fullSF = Upsample * obj.resUp; %Extent to which the ellipse images should be blown up.
            
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
            imgr = paintEllipse(imgr,obj.xCells,obj.yCells,obj.aCells/2,obj.lam*ones(size(obj.aCells))/2,-rad2deg(obj.thetCells),obj.cCells(:,1),1/fullSF,obj.boundConds,obj.xWidth,obj.yHeight);
            imgg = paintEllipse(imgg,obj.xCells,obj.yCells,obj.aCells/2,obj.lam*ones(size(obj.aCells))/2,-rad2deg(obj.thetCells),obj.cCells(:,2),1/fullSF,obj.boundConds,obj.xWidth,obj.yHeight);
            imgb = paintEllipse(imgb,obj.xCells,obj.yCells,obj.aCells/2,obj.lam*ones(size(obj.aCells))/2,-rad2deg(obj.thetCells),obj.cCells(:,3),1/fullSF,obj.boundConds,obj.xWidth,obj.yHeight);
            
            %Smooth to make nicer looking
            imgr = imgaussfilt(imgr,fullSF*obj.lam/4);
            imgg = imgaussfilt(imgg,fullSF*obj.lam/4);
            imgb = imgaussfilt(imgb,fullSF*obj.lam/4);
            
            outImg = cat(3,imgr,imgg,imgb);
%             imshow(outImg,'Parent',hand.CurrentAxes);
%             axis(hand.CurrentAxes,'tight')
            
            outImg = imresize(outImg,1/Downsample);
        end
    end
end