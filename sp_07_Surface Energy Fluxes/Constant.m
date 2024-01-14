%% --- Physical constants

physcon.vkc = 0.4;                              % von Karman constant
physcon.grav = 9.80665;                         % Gravitational acceleration (m/s2)
physcon.tfrz = 273.15;                          % Freezing point of water (K)
physcon.sigma = 5.67e-08;                       % Stefan-Boltzmann constant (W/m2/K4)
physcon.mmdry = 28.97 / 1000;                   % Molecular mass of dry air (kg/mol)
physcon.mmh2o = 18.02 / 1000;                   % Molecular mass of water (kg/mol)
physcon.cpd = 1005.0;                           % Specific heat of dry air at constant pressure (J/kg/K)
physcon.cpw = 1846.0;                           % Specific heat of water vapor at constant pressure (J/kg/K)
physcon.rgas = 8.31446;                         % Universal gas constant (J/K/mol)
physcon.cwat = 4188.0;                          % Specific heat of water (J/kg/K)
physcon.cice = 2117.27;                         % Specific heat ice (J/kg/K)
physcon.rhowat = 1000.0;                        % Density of water (kg/m3)
physcon.rhoice = 917.0;                         % Density of ice (kg/m3)
physcon.cvwat = physcon.cwat * physcon.rhowat;  % Heat capacity of water (J/m3/K)
physcon.cvice = physcon.cice * physcon.rhoice;  % Heat capacity of ice (J/m3/K)
physcon.tkwat = 0.57;                           % Thermal conductivity of water (W/m/K)
physcon.tkice = 2.29;                           % Thermal conductivity of ice (W/m/K)
physcon.hfus = 0.3337e6;                        % Heat of fusion for water at 0 C (J/kg)
physcon.hvap = 2.501e6;                         % Latent heat of evaporation (J/kg)
physcon.hsub = physcon.hfus + physcon.hvap;     % Latent heat of sublimation (J/kg)



solcon = 1364.0;          % Solar constant (W/m2)
