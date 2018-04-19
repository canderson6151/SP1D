%
% CreatePotentialOp
%
% This function creates the sparse matrix representation of an 
% operator that is the multiplication by a potential. 
%
% Inputs
%
% potential: A vector of length pointCount containing the values
%            of the potential function. The values at the first 
%            and last index are ignored.
function [P] = CreatePotentialOp(potential)
 

pointCount = length(potential);
M          = length(potential)-2;

% Create diagonal sparse matrix for the potential

iIndex = 1:M;
jIndex = 1:M;
P      = sparse(iIndex,jIndex,potential(2:pointCount-1),M,M);


 
