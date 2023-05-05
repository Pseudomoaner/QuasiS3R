function [fT,fR,fPar] = calcFrictionTensors(a,u,f0,fricAnis)
%CALCFRICTIONTENSORS calculates the translational and rotational friction 
%tensors for each rod.
%
%   INPUTS:
%       -a: The aspect ratios of all the rods
%       -u: The orientation unit vectors of all the rods
%       -f0: The Stokesian friction coefficient
%       -fricAnis: Specifies the extent of the anisotropy of the friction
%       tensor. Set to zero for isotropic, 1 for the original value as
%       used in the Wensink framework. Other values are inter/extrapolated
%       between these points.
%
%   OUTPUTS:
%       -fT: The translational friction tensor for each rod
%       -fR: The rotational friction tensor for each rod
%       -fPar: The parallel component of the friction tensor for each rod
%       (used in velocity calculations)
%
%   Author: Oliver J. Meacock, 2020

[fPar,fPerp,fRot] = calcGeomFactors(a);

meanF = (fPar + fPerp)/2; %Take the mean of the two components for each rod
fPerp = meanF + ((fPerp-fPar)*fricAnis)/2;
fPar = meanF - ((fPerp-fPar)*fricAnis)/2;

fR = f0 * fRot;
fT = zeros(3,3,size(a,1));
for i = 1:size(a,1)
    uSq = u(i,:)'*u(i,:);
    fT(:,:,i) = f0 * (fPar(i)*uSq + fPerp(i)*(eye(3) - uSq)); %Sets the basis of the different parallel and perpendicular friction coefficients to be oriented along the cell axis.
end