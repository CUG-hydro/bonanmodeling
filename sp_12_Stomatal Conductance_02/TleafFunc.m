function [flux, tleaf_dif] = TleafFunc (tleaf_val, varargin)

% Calculate leaf temperature and fluxes for an input leaf temperature
% (tleaf_val) and compare the new temperature to the prior
% temperature. This function returns a value tleaf_dif = 0 when leaf
% temperature does not change between iterations.

if length(varargin) == 1 && iscell(varargin{1}); varargin = varargin{1}; end
[physcon, atmos, leaf, flux] = flatten(varargin);

if (tleaf_val < 0)
   error ('TleafFunc error')
end

% --- Current value for leaf temperature
flux.tleaf = tleaf_val;

% --- Leaf boundary layer conductances
[flux] = LeafBoundaryLayer (physcon, atmos, leaf, flux);

% --- Leaf photosynthesis and stomatal conductance
[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);

% --- Leaf temperature and energy fluxes
[flux] = LeafTemperature (physcon, atmos, leaf, flux);

% --- Compare with prior value for leaf temperature
tleaf_dif = flux.tleaf - tleaf_val;
