% Supplemental program 7.1
import module_Radiation.*
import module_HydroTools.*
import module_Ipaper.*

solcon = 1364;                      % Solar constant (W/m2)

% --- Model run control parameters
nday = 30;                          % Number of days to simulate, repeating the same diurnal cycle
soilvar.method = 'excess-heat';     % Phase change: use 'excess-heat' or 'apparent-heat-capacity'

fluxvar.profiles = 'MOST';        % Use Monin-Obukhov similarity theory
% fluxvar.profiles = 'RSL';         % Use canopy coupling with roughness sublayer theory

fluxvar.bucket = 'no_bucket';     % Soil wetness factor = 1
% fluxvar.bucket = 'use_bucket';    % Use bucket model hydrology for soil wetness factor

% --- Atmospheric forcing at a reference height
forcvar.zref = 30.0;      % Reference height (m)
tmean = 25.0;             % Mean daily air temperature (C)
trange = 10.0;            % Temperature range for diurnal cycle (C)
relhum = 70.0;            % Relative humidity (%)
RH = 70;
forcvar.uref = 3.0;       % Wind speed at reference height (m/s)
forcvar.pref = 101325;    % Atmospheric pressure (Pa)
forcvar.rain = 0.0;       % Rainfall (kg H2O/m2/s)
forcvar.snow = 0.0;       % Snowfall (kg H2O/m2/s)
% [Ta, DRT, RH, U, P, RAIN, SNOW] = textread('forcing.txt','%f %f %f %f %f %f %f');

doy = 182.0;              % Day of year (1 to 365) for solar radiation
lat = 40.0 * pi/180;      % Latitude (degrees -> radians) for solar radiation

% --- Site characteristics
vis = 1; nir = 2;               % Waveband indices for visible and near-infrared
alb_surf(vis) = 0.10;           % Snow-free surface albedo for visible waveband (-)
alb_surf(nir) = 0.20;           % Snow-free surface albedo for near-infrared waveband (-)
alb_snow(vis) = 0.95;           % Snow albedo for visible waveband (-)
alb_snow(nir) = 0.70;           % Snow albedo for near-infrared waveband (-)
surfvar.emiss = 0.98;           % Surface emissivity (dimensionless)
snow_mask = 100.0;              % Snow albedo masking depth (kg H2O/m2)

soilvar.soil_texture = 5;       % Soil texture class

bucket.soil_water_max = 150.0;  % Maximum soil water (kg H2O/m2)
bucket.soil_beta_max = 0.75;    % Soil water at which soil_beta = 1 (fraction of soil_water_max)

surfvar.hc = 20.0;              % Canopy height (m)

%surfvar.hc = 0.5;
%forcvar.zref = 10.5;
surfvar.LAI = 5.0;              % Leaf area index (m2/m2)

% For Harman and Finnigan (2007, 2008) roughness sublayer (RSL)
lad = surfvar.LAI / surfvar.hc; % Leaf area density (m2/m3)
cd = 0.20;                      % Leaf drag coefficient (dimensionless)
surfvar.Lc = 1 / (cd * lad);    % Canopy density length scale (m)
surfvar.rc = 0.2;               % Leaf Nusselt number (heat) or Stanton number (scalar)

% For Monin-Obukhov parameterization (MOST)
fluxvar.disp = 0.67 * surfvar.hc;   % Displacement height (m)
fluxvar.z0m = 0.13 * surfvar.hc;    % Roughness length for momentum (m)
fluxvar.z0c = 0.10 * fluxvar.z0m;   % Roughness length for scalars (m)

%% --- Soil variables
% Number of layers in soil profile

% Soil layer thickness (m)
depths = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058];
soilvar = soil_depth_init(soilvar, depths);

% --- Initial conditions
% Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)
for i = 1:soilvar.nsoi
  % Temperature
  soilvar.tsoi(i) = tmean + physcon.tfrz;
  
  % Soil water at saturation (kg H2O/m2)
  h2osoi_sat = soilvar.watsat(soilvar.soil_texture) * physcon.rhowat * soilvar.dz(i);
  
  % Actual water content is some fraction of saturation. These are only used for soil
  % thermal properties and phase change. Note the inconsistency with the use of soil
  % water in the bucket model to calculate the soil wetness factor.
  satfrac = 0.85;
  if (soilvar.tsoi(i) > physcon.tfrz)
    soilvar.h2osoi_ice(i) = 0;
    soilvar.h2osoi_liq(i) = satfrac * h2osoi_sat;
  else
    soilvar.h2osoi_ice(i) = satfrac * h2osoi_sat;
    soilvar.h2osoi_liq(i) = 0;
  end
end

% Initial surface temperature (K) and vapor pressure (Pa)
fluxvar.tsrf = soilvar.tsoi(1);
[esat, desat] = satvap (fluxvar.tsrf - physcon.tfrz);
fluxvar.esrf = esat;

% Bucket model snow and soil water (kg H2O/m2)
bucket.snow_water = 0;
bucket.soil_water = bucket.soil_water_max;

%% --- Time stepping loop
% Main loop is NTIM iterations per day with a time step of DT seconds.
% This is repeated NDAY times to spinup the model from arbitrary initial
% soil temperatures.
dt = 1800;    % Time step (seconds)
ntim = round(86400/dt);