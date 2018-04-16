function [rho] = CreateDensity(energies,U,denParams)

% Extract parameters

layerCount     = denParams.layerCount;
pointCount     = denParams.pointCount;
E_fermi        = denParams.E_fermi;
hFactor        = denParams.hFactor;
effMassFactors = denParams.effMassFactors;
zLayerIndex    = denParams.zLayerIndex;

rho            = zeros(pointCount,1);

for k = 1:length(energies)
  if(energies(k) < E_fermi)
   for i = 1:layerCount
    iA = zLayerIndex(i);
    iB = zLayerIndex(i+1);
    rho(iA:iB) = (U(iA:iB,k).*U(iA:iB,k))*densityOfStates2D(energies(k),E_fermi,hFactor,effMassFactors(i));
   end

  % Overwrite values at internal layer boundaries with a value that is 
  % the average of the density of states on either side of the layer boundary.

   for i = 2:layerCount
    iStar = zLayerIndex(i);
    rho(iStar) =  0.5*( ...
      (U(iStar,k)*U(iStar,k))*densityOfStates2D(energies(k),E_fermi,hFactor,effMassFactors(i-1)) ...
    + (U(iStar,k)*U(iStar,k))*densityOfStates2D(energies(k),E_fermi,hFactor,effMassFactors(i))); 
   end
  end

  
end

