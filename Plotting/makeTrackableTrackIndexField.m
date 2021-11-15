function indexField = makeTrackableTrackIndexField(toMappings)
%MAKETRACKABLETRACKINDEXFIELD creates a new field for the trackableData
%strucure containing the index of each object, as indexed in toMappings.

indexField = cell(size(toMappings));

for i = 1:size(toMappings,2)
    indexField{i} = toMappings{i}(:,1);
end