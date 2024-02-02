function [flux] = StomataOptimization (physcon, atmos, leaf, flux)

% Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
% with water-use efficiency stomatal optimization

% --- Low and high initial estimates for gs (mol H2O/m2/s)

gs1 = 0.002;
gs2 = 2.0;

% --- Check for minimum stomatal conductance linked to low light or drought stress based
% on the water-use efficiency and cavitation checks for gs1 and gs2 (check1, check2)

[flux, check1] = StomataEfficiency (gs1, physcon, atmos, leaf, flux);
[flux, check2] = StomataEfficiency (gs2, physcon, atmos, leaf, flux);

if (check1 * check2 < 0)

   % Calculate gs using the function StomataEfficiency to iterate gs
   % to an accuracy of tol (mol H2O/m2/s)

   tol = 0.004;
   func_name = 'StomataEfficiency';
   [flux, root] = root_brent (func_name, gs1, gs2, tol, physcon, atmos, leaf, flux);
   flux.gs = root;

else

   % Low light or drought stress. Set gs to minimum conductance

   flux.gs = 0.002;

end

% --- Leaf fluxes for this gs

[flux] = LeafBoundaryLayer (physcon, atmos, leaf, flux);
[flux] = LeafTemperature (physcon, atmos, leaf, flux);
[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
[flux.psi_leaf] = LeafWaterPotential (physcon, leaf, flux);
