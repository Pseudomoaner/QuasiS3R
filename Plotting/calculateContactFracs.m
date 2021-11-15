function [stotContactFrac,ttosContactFrac] = calculateContactFracs(trackableData,fieldSettings)
%CALCULATECONTACTFRACS calculates the fraction of toxic to sensitive and
%sensitive to toxic contacts at each point in a simulation, relative to the
%total number of toxic and sensitive contacts respectively.

pxSize = 0.25;

%Construct a track index field for the overlap field painting function
for i = 1:size(trackableData.Force,2)
    trackableData.TrackIndex{i} = 1:size(trackableData.Force{i},1);
end
overlaps = paintOverlapFieldTracks(trackableData,pxSize,fieldSettings.fieldWidth,fieldSettings.fieldHeight);

stotContactFrac = zeros(size(trackableData.Force)); %Fraction of sensitive to toxic contacts to all sensitive contacts
ttosContactFrac = zeros(size(trackableData.Force)); %Fraction of toxic to sensitive contacts to all toxic contacts

for i = 1:size(trackableData.Force,2) %Loop over time
    fprintf('Frame %d of %d.\n',i,size(trackableData.Force,2))
    
    currOverlap = overlaps(:,:,i);

    stosCount = 0;
    stotCount = 0;
    ttotCount = 0;
    ttosCount = 0;
    for j = 1:size(trackableData.Force{i},1)
        currCell = makeExpandedRodProfileTrackableData(trackableData,j,i,fieldSettings.fieldWidth,fieldSettings.fieldHeight,pxSize);
        tgtList = unique(currOverlap(currCell(:)));
        tgtList(tgtList == 0) = [];
        tgtList(tgtList == j) = [];
        popList = trackableData.Population{i}(tgtList);

        if trackableData.Population{i}(j) == 's'
            stosCount = stosCount + sum(popList == 's');
            stotCount = stotCount + sum(popList == 't');
        else
            ttosCount = ttosCount + sum(popList == 's');
            ttotCount = ttotCount + sum(popList == 't');
        end
    end

    stotContactFrac(i) = stotCount/(stotCount + stosCount);
    ttosContactFrac(i) = ttosCount/(ttosCount + ttotCount);
end