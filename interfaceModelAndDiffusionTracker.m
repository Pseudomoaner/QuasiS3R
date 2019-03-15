function TrackableData = interfaceModelAndDiffusionTracker(Field,TrackableData,FrameNo)
%Creates an instance of PCs that is based on the given model frame.

if isempty(TrackableData)
    TrackableData = struct();
end

TrackableData.Length{FrameNo} = Field.aCells;
oriTmp = rad2deg(Field.thetCells);
oriTmp(oriTmp < -90) = oriTmp(oriTmp < -90) + 180;
oriTmp(oriTmp > 90) = oriTmp(oriTmp >90) + 180;
TrackableData.Orientation{FrameNo} = oriTmp;
TrackableData.Centroid{FrameNo} = [Field.xCells,Field.yCells];
TrackableData.Tilt{FrameNo} = Field.phiCells;