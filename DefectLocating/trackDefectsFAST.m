function [posDefTrackStruct,negDefTrackStruct,trackSettings] = trackDefectsFAST(pDposition,nDposition,pDorientation,nDorientation,fS,pS,dt)

%Loop over timepoints
for t = 1:size(pDorientation,2)
    %Each row of pDposition{t} and nDposition{t} contains the (x,y) coordinates of a single defect
    trackableDataPos.Centroid{t} = pDposition{t};
    trackableDataNeg.Centroid{t} = nDposition{t};
    %Similarly for defect orientation:
    trackableDataPos.Orientation{t} = pDorientation{t}/2; %Need to adjust orientations so they vary between -90 and 90 degrees to interface with FAST (will remove transformation later)
    trackableDataNeg.Orientation{t} = (3*nDorientation{t})/2;
    %SpareFeat1 contains the defect charge information. Additional features MUST be called SpareFeat1, SpareFeat2 etc. in the trackable data structure to be read by the tracking module
    trackableDataPos.SpareFeat1{t} = ones(size(pDposition{t},1),1);
    trackableDataNeg.SpareFeat1{t} = -ones(size(nDposition{t},1),1);
end

trackSettings.SpareFeat1 = 0;
trackSettings.Centroid = 1;
trackSettings.Orientation = 1;
trackSettings.Velocity = 0;
trackSettings.Length = 0;
trackSettings.Area = 0;
trackSettings.Width = 0;
trackSettings.noChannels = 0;
trackSettings.availableMeans = [];
trackSettings.availableStds = [];
trackSettings.MeanInc = [];
trackSettings.StdInc = [];
trackSettings.SpareFeat2 = 0;
trackSettings.SpareFeat3 = 0;
trackSettings.SpareFeat4 = 0;

trackSettings.incProp = pS.incProp;
trackSettings.tgtDensity = pS.tgtDensity;
trackSettings.gapWidth = 1;
trackSettings.maxFrame = size(pDorientation,2);
trackSettings.minFrame = 1;
trackSettings.minTrackLen = pS.minTrackLength;
trackSettings.frameA = 1;
trackSettings.statsUse = 'Centroid';
trackSettings.pseudoTracks = false; %Variable set to true in extractFeatureEngine.m if the 'tracks' have come from a single frame.

trackSettings.dt = dt;
trackSettings.pixSize = pS.pixSize;
trackSettings.maxX = fS.xWidth;
trackSettings.maxY = fS.yHeight;
trackSettings.maxF = trackSettings.maxFrame;

debugSet = true; %Prevents modal locking of progress bars

[linkStatsPos,featMatsPos,featureStructPos,possIdxPos] = gatherLinkStats(trackableDataPos,trackSettings,debugSet);
[linkStatsNeg,featMatsNeg,featureStructNeg,possIdxNeg] = gatherLinkStats(trackableDataNeg,trackSettings,debugSet);

%Build feature matrices
[featMatsPos.lin,featMatsPos.circ] = buildFeatureMatricesRedux(trackableDataPos,featureStructPos,possIdxPos,trackSettings.minFrame,trackSettings.maxFrame);
[featMatsNeg.lin,featMatsNeg.circ] = buildFeatureMatricesRedux(trackableDataNeg,featureStructNeg,possIdxNeg,trackSettings.minFrame,trackSettings.maxFrame);

[TracksPos,InitialsPos] = doDirectLinkingRedux(featMatsPos.lin,featMatsPos.circ,featMatsPos.lin,featMatsPos.circ,linkStatsPos,trackSettings.gapWidth,false,debugSet);
[TracksNeg,InitialsNeg] = doDirectLinkingRedux(featMatsNeg.lin,featMatsNeg.circ,featMatsNeg.lin,featMatsNeg.circ,linkStatsNeg,trackSettings.gapWidth,false,debugSet);

trackDataNames = fieldnames(trackableDataPos);
rawTracksPos = struct();
for i = 1:size(trackDataNames,1)
    if i == 1
        [rawTracksPos.(trackDataNames{i}),trackTimesPos,rawToMappingsPos,rawFromMappingsPos] = extractDataTrack(TracksPos,InitialsPos,trackableDataPos.(trackDataNames{i})(trackSettings.minFrame:trackSettings.maxFrame),true);
        [rawTracksNeg.(trackDataNames{i}),trackTimesNeg,rawToMappingsNeg,rawFromMappingsNeg] = extractDataTrack(TracksNeg,InitialsNeg,trackableDataNeg.(trackDataNames{i})(trackSettings.minFrame:trackSettings.maxFrame),true);
    else
        rawTracksPos.(trackDataNames{i}) = extractDataTrack(TracksPos,InitialsPos,trackableDataPos.(trackDataNames{i})(trackSettings.minFrame:trackSettings.maxFrame),false);
        rawTracksNeg.(trackDataNames{i}) = extractDataTrack(TracksNeg,InitialsNeg,trackableDataNeg.(trackDataNames{i})(trackSettings.minFrame:trackSettings.maxFrame),false);
    end
end

[posDefTracks,fromPosMappings,toPosMappings] = processTracks(rawTracksPos,rawFromMappingsPos,rawToMappingsPos,trackSettings,trackTimesPos,debugSet);
[negDefTracks,fromNegMappings,toNegMappings] = processTracks(rawTracksNeg,rawFromMappingsNeg,rawToMappingsNeg,trackSettings,trackTimesNeg,debugSet);

%Assign tracks to two separate populations so you can easily use FAST's
%plotting methods
for i = 1:size(posDefTracks,2)
    posDefTracks(i).population = 1;
    posDefTracks(i).phi = posDefTracks(i).phi*2; %Reverse the transformation you applied in line 8
end
for i = 1:size(negDefTracks,2)
    negDefTracks(i).population = 2;
    negDefTracks(i).phi = negDefTracks(i).phi*2/3;
    %         procDefTracks(i).phi = procDefTracks(i).phi + (floor(rand(size(procDefTracks(i).phi))*3)-1)*120; %Randomise orientation in 120 degree increments at every timepoint
    negDefTracks(i).phi = wrapTo180(rad2deg(unwrap(deg2rad(negDefTracks(i).phi*3)))/3 + floor((rand(1)*3)-1)*120); %Randomise orientation in 120 degree increments at a single timepoint and then keep orientation fixed
end

posDefTrackStruct.tracks = posDefTracks;
posDefTrackStruct.toMappings = toPosMappings;
posDefTrackStruct.fromMappings = fromPosMappings;

negDefTrackStruct.tracks = negDefTracks;
negDefTrackStruct.toMappings = toNegMappings;
negDefTrackStruct.fromMappings = fromNegMappings;

debugprogressbar(1)