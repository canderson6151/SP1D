%
% Function SPoperator
%
% This function return phi_(n+1) that is computed using 
%
% phi_(n+1)  = -Delta^(-1) rho( phi_n, phiBase) - phi_n
%
% Where rho(phi_n, phiBase) is the density of electrons associated 
% with a potential given by phi = phi_n + phiBase
%
% Using Forward Euler with this operator with timestep dt is equivalent
% to using simple fixed point iteration with a relaxation parameter dt. 
%
%
%
%


function [phiNP1] = SPoperator(phiN,phiBase,potParams,denParams,L,K,D,sqrtInvD)

% Set Schroedinger operator potential

P  = CreatePotentialOp(phiN + phiBase);
S  = K + P;

% Create eigensystem

[eigVectors,eigVal] = CreateEigenSystem(S,sqrtInvD);

% Create density 

rho = CreateDensity(eigVal,eigVectors,denParams);

phiA    = 0.0;
phiB    = 0.0;
rhoStar = -rho;

phiNP1  = CreatePotential(rhoStar,phiA,phiB,potParams,L,D) - phiN;
