clear, clc;
import module_Radiation.*
import module_HydroTools.*
import module_Ipaper.*

sp_forcing
open forcvar
open surfvar
open soilvar
open fluxvar
%% 30天之后的稳态
j = 1;
i = 20;
% for j = 1:nday
fprintf('day = %6.0f\n',j)

% for i = 1:ntim
hour = i * (dt/86400 * 24); % Hour of day (0 to 24)

% Air temperature (K): use a sine wave with max (tmean + 1/2 trange) at 1400
% and min (tmean - 1/2 trange) at 0200
forcvar.tref = tmean + 0.5 * trange * sin(2*pi/24 * (hour-8)) + physcon.tfrz; % Ta_ref

% Vapor pressure (Pa) using constant relative humidity
[esat, desat] = satvap (forcvar.tref-physcon.tfrz);
forcvar.eref = (RH / 100) * esat;

[forcvar.cpair, forcvar.thref, forcvar.thvref, forcvar.mmair, forcvar.rhomol] = ...
  cal_Cp(forcvar.tref, forcvar.eref, forcvar.pref, forcvar.zref);

% Solar radiation at top of the atmosphere
[Rs_toa, Rs, Rs_dir, Rs_dif, coszen] = cal_Rs_toa(lat, doy, hour);
forcvar.solrad = Rs;                    % Total at surface

% Longwave radiation (W/m2)
forcvar.lwdown = (0.398e-05 * forcvar.tref^2.148) * physcon.sigma * forcvar.tref^4;

% Effective surface albedo is weighted combination of snow-free and
% snow albedos
fsno = bucket.snow_water / (bucket.snow_water + snow_mask);
alb_eff(vis) = alb_surf(vis) * (1 - fsno) + alb_snow(vis) * fsno;
alb_eff(nir) = alb_surf(nir) * (1 - fsno) + alb_snow(nir) * fsno;

% Radiative forcing: absorbed solar + incident longwave. This partitions
% solar radiation into 50% visible and 50% near-infrared wavebands.
fluxvar.qa = (1-alb_eff(vis)) * 0.5*forcvar.solrad ...
  + (1-alb_eff(nir)) * 0.5*forcvar.solrad + surfvar.emiss * forcvar.lwdown;

% Canopy conductance (mol/m2/s) - use a weighted average of sunlit and shaded leaves
gcan_min = 0.05;                           % Minimum conductance (mol/m2/s)
gcan_max = 0.2;                            % Maximum conductance (mol/m2/s)

ext = 0.5 / max(coszen, 0.0001);           % Light extinction coefficient
fsun = (1 - exp(-ext*surfvar.LAI)) / (ext*surfvar.LAI);  % Sunlit fraction of canopy

if (forcvar.solrad > 0)
  surfvar.gcan = (fsun * gcan_max + (1 - fsun) * gcan_min) * surfvar.LAI;
else
  surfvar.gcan = gcan_min * surfvar.LAI;
end

% Thermal conductivity and heat capacity
[soilvar] = soil_thermal_properties (physcon, soilvar);

% % Calculate the soil temperatures and surface fluxes
[fluxvar, soilvar, bucket] = surface_fluxes(physcon, forcvar, surfvar, soilvar, fluxvar, bucket, dt);

% [fluxvar, dx] = most(10, physcon, forcvar, surfvar, fluxvar)
% open fluxvar
