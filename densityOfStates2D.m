%
% A density of states function for a T = 0 (0 temperature). 
%
%

function density = densityOfStates2D(E,E_fermi,hFactor,effMassFactor)
%
% Formula used (double occupancy)
% 
% Density(E) = 2*((E_fermi -E)*effMass)/(2*pi*hbar^2)   E < E_fermi
%
% When expressed in terms of hFactor and effMassFactor, this becomes
%
% (1/(2*pi))*(effmassFactor/hFactor)
% 
density             = 0.0;
if(E < E_fermi) 
  density = abs(E_fermi-E)*((effMassFactor/hFactor)/(2.0*pi));
end
