%
% Function CreatePotential
%
% This function uses the sparse matrix representation 
% of a finite volume discretization of Poisson's problem 
% to solve the Poisson problem with Dirichlet boundary
% conditions. 
%
% Inputs
%
% f: right hand side of Poisson's equation. A vector of length
%    pointCount. The values at the first and last index corresponding
%    to values at boundary points are ignored.
%
% phiA, phiB : values specifying the boundary values of the solution
%              at the endpoints of the interval.
%
% potParams, L, D : Grid parameters and sparse matrix representations of the
%                   discrete approxiamtion to Poisson's equation. 
%
%                   The L and D sparse matrices are those created by 
%                   the CreateLapOp function. 
%



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
