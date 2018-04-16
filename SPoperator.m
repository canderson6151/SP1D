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
