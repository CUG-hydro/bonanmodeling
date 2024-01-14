classdef physcon
  %% --- Physical constants
  
  properties (Constant)
    vkc    = 0.4;                            % von Karman constant
    grav   = 9.80665;                        % Gravitational acceleration (m/s2)
    tfrz   = 273.15;                         % Freezing point of water (K)
    sigma  = 5.67e-08;                       % Stefan-Boltzmann constant (W/m2/K4)
    mmdry  = 28.97 / 1000;                   % Molecular mass of dry air (kg/mol)
    mmh2o  = 18.02 / 1000;                   % Molecular mass of water (kg/mol)
    cpd    = 1005.0;                         % Specific heat of dry air at constant pressure (J/kg/K)
    cpw    = 1846.0;                         % Specific heat of water vapor at constant pressure (J/kg/K)
    rgas   = 8.31446;                        % Universal gas constant (J/K/mol)
    cwat   = 4188.0;                         % Specific heat of water (J/kg/K)
    cice   = 2117.27;                        % Specific heat ice (J/kg/K)
    rhowat = 1000.0;                         % Density of water (kg/m3)
    rhoice = 917.0;                          % Density of ice (kg/m3)
    cvwat  = physcon.cwat * physcon.rhowat;  % Heat capacity of water (J/m3/K)
    cvice  = physcon.cice * physcon.rhoice;  % Heat capacity of ice (J/m3/K)
    tkwat  = 0.57;                           % Thermal conductivity of water (W/m/K)
    tkice  = 2.29;                           % Thermal conductivity of ice (W/m/K)
    hfus   = 0.3337e6;                       % Heat of fusion for water at 0 C (J/kg)
    hvap   = 2.501e6;                        % Latent heat of evaporation (J/kg)
    hsub   = physcon.hfus + physcon.hvap;    % Latent heat of sublimation (J/kg)
    solcon = 1367; % W/m2
  end

end
