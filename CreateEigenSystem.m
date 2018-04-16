%
% Function: CreateEigenSystem
%
% This m-file creates computes all the eigenvectors and 
% eigenvectors associated with the generalied eigenproblem
% 
%    S u = lambda*D u
%
% and D is a diagonal matrix with positive entries.
%
% Procedure: Convert the problem into a standard eigenproblem
% and then use the eig(...) routine. 
%
% Letting v = D^(1/2) u the eigenproblem becomes 
%
% D^(-1/2) * S * D^(-1/2) v = lambda v 
%
% then u = D^(-1/2) v are the eigenvectors 
% and orthonomality of the eigenvectors v in the standard 
% inner product leads to orthogonality of the u vectors in the
% mesh scaled inner product, e.g. <u_j, D u_i> = delta_(i,j) 
%
% The eigenvectors returned by this method are defined at all grid 
% points (the endpoint values are set to 0). The eigenvalues and
% corresponding eigenvectors are sorted in ascending order. 
%
function [eigVectors,eigVal] = CreateEigenSystem(S,sqrtInvD)
%
% Convert to a standard symmetric eigenvalue problem 
% 
A                      = sqrtInvD*S*sqrtInvD;
[eigVecTemp,eigValues] = eig(full(A));


% Scale appropriately. Note, the resulting eigVectors will NOT be 
% orthonormal with respect to the standard inner product, however, 
% they will be orthonormal with respect to the mesh weighted innner product 
% (which is what one wants). 
%

eigVecTemp = sqrtInvD*eigVecTemp;

% Create sorting index 

[d,ind] = sort(diag(eigValues),'ascend');

eigValTemp   = diag(eigValues);
eigVal       = eigValTemp(ind); 

M          = size(S,1);
pointCount = M + 2;

% Add in boundary values to sorted eigenvectors 

eigVectors  = zeros(pointCount,M);
eigVectors(2:pointCount-1,:) = eigVecTemp(:,ind);

