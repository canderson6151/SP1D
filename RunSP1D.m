%
% A sample program written using Matlab syntax that
% determines the self-consistent solution of the 
% 1D Schroedinger-Poisson equation. 
%
% Assumptions:
%
% The boundary conditions for the potential are specified
%
% Homogeneous Dirichlet boundary conditions for the eigenfunctions 
% of the Schroedinger operator 
%
% Units : distance in nm,
%       : energy   in eV
%       : time     in sec
% 
% The coefficient of the kinetic energy operator 
% is (hbar^2/(2*mStar)). where Mstar = effMassFactor * electron rest mass.
% In the units used, this value is specified by the relation
% (hbar^2/(2*mStar)) = hFactor/(effMassFactor) with hFactor = 0.0381
%
% Created for Field Institute Mini-Course
% "Numerical Techniques for Simulating Collective and Coherent Qusntum States"
%
% Author  : Chris Anderson
% Date    : April 15, 2018
% License : GPLv3 (you can do what you want with it)
%
% An "alpha" version (necessary validation still needs to be done)
%
% Note: Since Dirichlet conditions are used, the size of the 
% sparse matrix representation of the operators is 
% pointCount - 2 X pointCount - 2 where pointCount is the total number 
% of grid points. However, for the functions that use these 
% sparse matrices and have vector inputs (for right hand sides, etc.) 
% the lengths of the input vectors is pointCount.  
%

E_fermi        = 0.0;
 
phiA           = 0.0;
phiB           = 0.0;

Mat_1_dC  = 14.0; % dielectric constant
Mat_1_eM  = 0.06; % effective mass factor 
Mat_1_bS  = -0.2; % band shift

Mat_2_dC =  14.0; % dielectric constant
Mat_2_eM  = 0.07; % effective mass factor 
Mat_2_bS  = 0.0;  % band shift

hFactor                    =  0.0381;
dielectricConversionFactor = .055271;

diElecConstants  = [Mat_2_dC, Mat_1_dC, Mat_2_dC, Mat_1_dC, Mat_2_dC];
effMassFactors   = [Mat_2_eM, Mat_1_eM, Mat_2_eM, Mat_1_eM, Mat_2_eM];
bandShift        = [Mat_2_bS, Mat_1_bS, Mat_2_bS, Mat_1_bS, Mat_2_bS];  

diElecCoeff      =  dielectricConversionFactor*diElecConstants;
kineticCoeff     = -hFactor./effMassFactors;  
           
%          
% Set up layer geometry and set up grid (widths currently in nm)
%   

layerWidths     = [30.0, 10.0, 10.0, 12.0, 30.0];
panelCounts     = 4*[ 10 ,  10,  10,   10,    10]; 

pointCount      = sum(panelCounts)+1;

layerCount      = length(panelCounts);
layerMeshSizes  = layerWidths./panelCounts;
zLayerIndex     = zeros(1,layerCount+1);

% Set up indices of grid points corresponding to layer boundaries.

zLayerIndex  = zeros(1,layerCount+1);
zLayerIndex(1) = 1;
for i = 1:layerCount
  zLayerIndex(i+1) = zLayerIndex(i) + panelCounts(i);
end 
 
% Set up z-coordinate locations of grid points 

zGrid         = zeros(pointCount,1);  
zIndex        = 1;
zGrid(zIndex) = 0;

for i = 1:layerCount
  hz = layerMeshSizes(i);
  for j = 2:(panelCounts(i)+1)
    zIndex = zIndex+1;
    zGrid(zIndex) = zGrid(zIndex-1) + hz;
   end
end
    
%
% Create band-shift potential. This potential is discontinuous across the 
% layer boundaries, so to keep it simple, we set the value on the layer 
% boundary to be the average of the values on either side of the boundary. 
%

bandShiftPotential   = zeros(pointCount,1);  

for i = 1:layerCount
  bandShiftPotential(zLayerIndex(i):zLayerIndex(i+1)) = bandShift(i);
end

% Fix up at internal layer boundaries by specifying a value that is 
% the average of the band shift on either side of the layer boundary.

for i = 2:layerCount
  bandShiftPotential(zLayerIndex(i)) = (bandShift(i) + bandShift(i-1))/2.0;
end

% Create potential operator components 

potParams.pointCount     = pointCount;
potParams.layerCount     = layerCount;
potParams.coeff          = diElecCoeff;
potParams.zLayerIndex    = zLayerIndex;
potParams.layerMeshSizes = layerMeshSizes;

[L, D, sqrtD, sqrtInvD] = CreateLapOp(potParams);

% Create kinetic operator component  

kineticParams.pointCount     = pointCount;
kineticParams.layerCount     = layerCount;
kineticParams.coeff          = kineticCoeff;
kineticParams.zLayerIndex    = zLayerIndex;
kineticParams.layerMeshSizes = layerMeshSizes;

[K, D, sqrtD, sqrtInvD] = CreateLapOp(kineticParams);

%
% Density operator parameters 
%

denParams.layerCount     = layerCount;
denParams.pointCount     = pointCount;
denParams.E_fermi        = E_fermi;
denParams.hFactor        = hFactor;
denParams.effMassFactors = effMassFactors;
denParams.zLayerIndex    = zLayerIndex;


% Capture weights for mesh weighted inner products for grid functions
% defined at all grid points. 

meshWeight      = zeros(1,pointCount);
meshWeight(1)   = layerMeshSizes(1)/2.0;
for i = 2:pointCount-1
  meshWeight(i) = D(i-1,i-1);
end
meshWeight(pointCount) = layerMeshSizes(layerCount)/2.0;

%
% Carry out 10 self-consistent iterations with relaxation
% factor dt. This self-consistent iteration is implemented
% as the Forward Euler method applied to the equation 
%
% dphiTilde/dt = -DELTA^(-1) rho_e( phiTilde, phiBase) - phiTilde
%
% See Notes Lecture 4
%
% A fixed number of timesteps is taken for demonstration,
% typically a stopping condition dependent upon the 
% difference between iterates would be used (along with
% a condition that stops if too many iterations have 
% been taken)
%

rhoZero   = zeros(pointCount,1);
phiBase   = CreatePotential(rhoZero,phiA,phiB,potParams,L,D);
phiBase   = phiBase + bandShiftPotential;

phiN = zeros(pointCount,1);

dt   = 0.1;
diff = zeros(pointCount,1);
for k = 1:10
  phiNP1 = phiN + dt*SPoperator(phiN,phiBase,potParams,denParams,L,K,D,sqrtInvD);
  norm((phiNP1-phiN)/dt,inf)
  phiN   = phiNP1;
end

phiFinal = phiN + phiBase;

%
% Create eigensystem and density associated with self-consistent 
% solution.
%

P  = CreatePotentialOp(phiFinal);
S  = K + P;
[eigVectors,eigVal] = CreateEigenSystem(S,sqrtInvD);
rho = CreateDensity(eigVal,eigVectors,denParams);

% Plot the final potential 

plot(zGrid,phiFinal);

pause

% Plot the density 

plot(zGrid,rho);







