function [phi] = CreatePotential(f,phiA,phiB,potParams,L,D)

% Extract parameters 

pointCount = potParams.pointCount;
layerCount = potParams.layerCount;
dCoeffA    = potParams.coeff(1);
dCoeffB    = potParams.coeff(layerCount);
hzA        = potParams.layerMeshSizes(1);
hzB        = potParams.layerMeshSizes(layerCount);

M          = pointCount-2;
% 
% Copy over right hand side 
%
fStar      = f(2:pointCount-1);

% Add boundary forcing to right hand side

fStar(1) = fStar(1) - (phiA*dCoeffA)/hzA;
fStar(M) = fStar(M) - (phiB*dCoeffB)/hzB;

% Solve linear system 

phiStar    = L\(D*fStar);

% Add in boundary values

phi    = zeros(pointCount,1);

phi(1)              = phiA;
phi(2:pointCount-1) = phiStar;
phi(pointCount)     = phiB;
