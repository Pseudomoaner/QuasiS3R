function [rmsd,times] = applyRmsd(maxF,data,dt)
%APPLYRMSD applies the RMSD analysis pipeline to the given input data.
%Assumes format is the same as FAST/TrackMate.
%
%   INPUTS:
%       -maxF: The maximal frame index
%       -data: The input track data
%       -dt: The size (in pysical units) of the timestep.
%
%   OUTPUTS:
%       -rmsd: The RMSD function for the x and y coordinates of the input
%       data
%       -times: The times (in physical units) of the RMSD samples.
%
%Author: Oliver J. Meacock, (c) 2020

sums = zeros(maxF,1);
counts = zeros(maxF,1);

for j = 1:size(data,2)
    for k = 1:length(data(j).x) %Treat each sub-track as a separate track
        currCents = [data(j).x(k:end),data(j).y(k:end)];
        currTimes = data(j).times(k:end); %Time indices, rather than real (minute) values.
        currTimes = currTimes - currTimes(1) + 1; %Center time points on start of track.
        
        currXDisps = currCents(:,1) - currCents(1,1);
        currYDisps = currCents(:,2) - currCents(1,2);
        
        sums(currTimes) = sums(currTimes) + currXDisps.^2 + currYDisps.^2;
        counts(currTimes) = counts(currTimes) + 1;
    end
end

rmsd = sqrt(sums./counts);
times = linspace(0,maxF*dt,maxF);