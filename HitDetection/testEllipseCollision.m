function collision = testEllipseCollision(inField,i,j)
%TESTELLIPSECOLLISION tests whether two input ellipses are colliding using
%a 1D optimisation approach. From
%https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect.
%
%   INPUTS:
%       -inField: An instance of a WensinkField object
%       -i,j: The indices of the rods you would like to check the overlap
%       for.
%
%   OUTPUTS:
%       -collision: True if ellipses are colliding, false if not.

%Get unit orientation vectors
ux = cos(inField.thetCells([i,j]));
uy = sin(inField.thetCells([i,j]));

%Deal with periodic BCs if needed
if strcmp(inField.boundConds,'periodic')
    x = inField.xCells([i,j]);
    if abs(x(1) - x(2)) > inField.xWidth/2
        if x(1) < x(2)
            x(1) = x(1) + inField.xWidth;
        else
            x(2) = x(2) + inField.xWidth;
        end
    end
    y = inField.yCells([i,j]);
    if abs(y(1) - y(2)) > inField.yHeight/2
        if y(1) < y(2)
            y(1) = y(1) + inField.yHeight;
        else
            y(2) = y(2) + inField.yHeight;
        end
    end
else
    x = inField.xCells([i,j]);
    y = inField.yCells([i,j]);
end

%Calculate effective aspect ratios
maj = inField.lam*inField.aCells([i,j])/2 + [inField.contactRange,0];
min = [inField.lam/2,inField.lam/2] + [inField.contactRange,0];

a = [x(1);y(1)];
b = [x(2);y(2)];
A = [ux(1),uy(1);uy(1),-ux(1)] * [1/(maj(1)^2), 0; 0, 1/(min(1)^2)] * [ux(1),uy(1);uy(1),-ux(1)];
B = [ux(2),uy(2);uy(2),-ux(2)] * [1/(maj(2)^2), 0; 0, 1/(min(2)^2)] * [ux(2),uy(2);uy(2),-ux(2)];

Ainv = inv(A);
Binv = inv(B);
centDiff = b-a;

testFunc = @(s)1-centDiff'*inv(((1/(1-s))*Ainv + (1/s)*Binv))*centDiff;

[~,fval] = fminbnd(testFunc,0,1);

collision = fval >= 0;