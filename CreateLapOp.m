% Function CreateLapOp
%
% Set up the symmetric tri-diagonal system of equations obtained from a 
% finite volume approximations to the Laplace operator with 
% Dirichlet boundary values. The underlying grid is a semi-layered grid --
% a grid that consists of a collection of layers, each layer of which
% has a uniform grid. 
%
% The diffusion coefficient is assumed to be constant over each
% layer, but different layers may have differnt constants. 
%
% Since the boundary values are given, the linear system is 
% a (pointCount-2) X (pointCount-2) tri-diagonal matrix.
%
% The matrices constructed by this routine are used by the
% CreatePotential m-file, as well as by the m-files used to
% compute the eigensystem of the Schroedinger operator. 
%
function [A,D,sqrtD,sqrtInvD] = CreateLapOp(potParams)

M              = potParams.pointCount - 2;
zLayerIndex    = potParams.zLayerIndex;
layerMeshSizes = potParams.layerMeshSizes;
coeff          = potParams.coeff;
layerCount     = length(layerMeshSizes);

A        = sparse(M,M);
D        = sparse(M,M);
sqrtD    = sparse(M,M);
sqrtInvD = sparse(M,M);

% Create row index bounds for rows associated with equations at points internal
% to the layers.

matLayerIndex(2:layerCount)   = zLayerIndex(2:layerCount) - ones(1,layerCount-1);
matLayerIndex(1)             = 1;
matLayerIndex(layerCount+1)  = M;

% First equation adjacent 

hz           =   layerMeshSizes(1);
i            =   matLayerIndex(1);

A(i,i)   = -(2.0*coeff(1))/(hz);
A(i,i+1) =      (coeff(1))/(hz);
D(i,i)        = hz;
sqrtD(i,i)    = sqrt(hz);
sqrtInvD(i,i) = 1.0/sqrt(hz);

for k = 1:layerCount
    hz           =   layerMeshSizes(k);
    for i = matLayerIndex(k)+ 1 : matLayerIndex(k+1)-1
      A(i,i-1) =      (coeff(k))/(hz);
      A(i,i)   = -(2.0*coeff(k))/(hz);
      A(i,i+1) =      (coeff(k))/(hz);
      
      D(i,i)        = hz;
      sqrtD(i,i)    = sqrt(hz);
      sqrtInvD(i,i) = 1.0/sqrt(hz);
    end
    
    % equations at internal layer boundary points 
    
    if(k < layerCount) 
      hzA          = layerMeshSizes(k);
      hzB          = layerMeshSizes(k+1);
      i            =  matLayerIndex(k+1);
      A(i,i-1)     =    coeff(k)/(hzA);
      A(i,i)       = -( coeff(k)/(hzA) + coeff(k+1)/(hzB)); 
      A(i,i+1)     =    coeff(k+1)/(hzB); 
      
      D(i,i)        = (hzA + hzB)/2.0;
      sqrtD(i,i)    = sqrt((hzA + hzB)/2.0);
      sqrtInvD(i,i) = 1.0/sqrt((hzA + hzB)/2.0);
    end
end


% Last equation   

k            =   layerCount;
hz           =   layerMeshSizes(k);
i            =   matLayerIndex(k+1);

A(i,i-1) =      (coeff(k))/(hz);
A(i,i)   = -(2.0*coeff(k))/(hz);
A(i,i-1) =      (coeff(k))/(hz);

D(i,i)        = hz;
sqrtD(i,i)    = sqrt(hz);
sqrtInvD(i,i) = 1.0/sqrt(hz);
