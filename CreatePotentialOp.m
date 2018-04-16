

function [P] = CreatePotentialOp(potential)
 

pointCount = length(potential);
M          = length(potential)-2;

% Create diagonal sparse matrix for the potential

iIndex = 1:M;
jIndex = 1:M;
P      = sparse(iIndex,jIndex,potential(2:pointCount-1),M,M);


 