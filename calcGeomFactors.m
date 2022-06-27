function [fPar,fPerp,fRot] = calcGeomFactors(a)
%CALCGEOMFACTORS calculates the geometrical contributions to the friction
%components.
%
%   INPUTS:
%       -a: The aspect ratio of a given cell.
%
%   OUTPUTS:
%       -fPar: The component of the friction tensor parallel to the rod's
%       long axis
%       -fPerp: The component of the friction tensor perpendicular to the
%       rod's long axis
%       -fRot: The rotational friction coefficient
%
%   Author: Oliver J. Meacock, 2020

fPar = (2 * pi)./(log(a) - 0.207 + (0.980./a) - (0.133./(a.^2)));
fPerp = (4 * pi)./(log(a) + 0.839 + (0.185./a) + (0.233./(a.^2)));
fRot = ((a.^2) * pi)./(3*(log(a) - 0.662 + (0.917./a) - (0.050./(a.^2))));
