function soilvar = soil_depth_init(soilvar, dz)
% %% INPUTS
% - dz              : depth of each layer
% 
% %% OUTPUTS
% - z               : z_i
% - z_plus_onehalf  : z_{i+1/2}
% - dz_plus_onehalf : dz_{i+1/2} = z_{i} - z_{i+1}

% dz = ones(1, nsoi) * 0.025; % depth of each lasyer (m)


%% Figure 5.3, soil_depth_init
% Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
% z_{i+1/2}
nsoi = length(dz);

z = zeros(1, nsoi);
z_plus_onehalf = zeros(1, nsoi);
dz_plus_onehalf = zeros(1, nsoi);

z_plus_onehalf(1) = -dz(1);
for i = 2:nsoi
  z_plus_onehalf(i) = z_plus_onehalf(i-1) - dz(i);
end

% Soil depth (m) at center of layer i (negative distance from surface)
z(1) = 0.5 * z_plus_onehalf(1);
for i = 2:nsoi
  z(i) = 0.5 * (z_plus_onehalf(i-1) + z_plus_onehalf(i));
end

% Thickness between between z(i) and z(i+1)
for i = 1:nsoi-1
  dz_plus_onehalf(i) = z(i) - z(i+1);
end
dz_plus_onehalf(nsoi) = 0.5 * dz(nsoi);

% OUTPUT
soilvar.nsoi = nsoi;

soilvar.dz              = dz;
soilvar.z               = z;
soilvar.z_plus_onehalf  = z_plus_onehalf;
soilvar.dz_plus_onehalf = dz_plus_onehalf;

%% --- Initialize soil texture variables
% Soil texture classes (Cosby et al. 1984. Water Resources Research 20:682-690)
%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silty loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay
soilvar.silt = [ 5.0, 12.0, 32.0, 70.0, 39.0, 15.0, 56.0, 34.0,  6.0, 47.0, 20.0]; % Percent silt
soilvar.sand = [92.0, 82.0, 58.0, 17.0, 43.0, 58.0, 10.0, 32.0, 52.0,  6.0, 22.0]; % Percent sand
soilvar.clay = [ 3.0,  6.0, 10.0, 13.0, 18.0, 27.0, 34.0, 34.0, 42.0, 47.0, 58.0]; % Percent clay

% Volumetric soil water content at saturation (porosity)
% (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
soilvar.watsat = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482];

end
