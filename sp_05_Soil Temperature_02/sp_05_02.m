%% Supplemental program 5.2
% Diurnal cycle of soil temperature with phase change using
% "excess heat" or "apparent heat capacity"

%% --- Physical constants in physcon structure
physcon.tfrz = 273.15;                         % Freezing point of water (K)
physcon.cwat = 4188.0;                         % Specific heat of water (J/kg/K)
physcon.cice = 2117.27;                        % Specific heat of ice (J/kg/K)
physcon.rhowat = 1000.0;                       % Density of water (kg/m3)
physcon.rhoice = 917.0;                        % Density of ice (kg/m3)
physcon.cvwat = physcon.cwat * physcon.rhowat; % Heat capacity of water (J/m3/K)
physcon.cvice = physcon.cice * physcon.rhoice; % Heat capacity of ice (J/m3/K)
physcon.tkwat = 0.57;                          % Thermal conductivity of water (W/m/K)
physcon.tkice = 2.29;                          % Thermal conductivity of ice (W/m/K)
physcon.hfus = 0.3337e6;                       % Heat of fusion for water at 0 C (J/kg)

%% --- Model run control parameters
tmean = physcon.tfrz + 15.0;       % Mean daily air temperature for diurnal cycle (K)
trange = 10.0;                     % Temperature range for diurnal cycle (K)
nday = 200;                        % Number of days
soilvar.soil_texture = 1;          % Soil texture class: sand
soilvar.method = 'excess-heat';    % Use excess heat for phase change

%soilvar.soil_texture = 11;         % Soil texture class: clay
%soilvar.method = 'apparent-heat-capacity'; % Use apparent heat capacity for phase change

% --- Initialize soil layer variables
% Soil layer thickness (m)
dz = ones(1, 120) * 0.025; % depth of each lasyer (m), 120 layers
soilvar = soil_depth_init(soilvar, dz);

% Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)
for i = 1:soilvar.nsoi
  % Temperature (K)
  soilvar.tsoi(i) = physcon.tfrz + 2.0;
  
  % Soil water at saturation (kg H2O/m2)
  h2osoi_sat = soilvar.watsat(soilvar.soil_texture) * physcon.rhowat * soilvar.dz(i);
  
  % Actual water content is some fraction of saturation
  if (soilvar.tsoi(i) > physcon.tfrz)
    soilvar.h2osoi_ice(i) = 0;
    soilvar.h2osoi_liq(i) = 0.8 * h2osoi_sat;
  else
    soilvar.h2osoi_ice(i) = 0.8 * h2osoi_sat;
    soilvar.h2osoi_liq(i) = 0;
  end
end

% --- Time stepping loop to increment soil temperature
% Main loop is NTIM iterations per day with a time step of DT seconds.
% This is repeated NDAY times.

dt = 1800;      % Time step (seconds), 0.5hour
ntim = round(86400/dt);

m = 0; % Counter for output file
for iday = 1:nday
  fprintf('day = %6.0f\n',iday)
  
  for itim = 1:ntim
    % Hour of day
    hour = itim * (dt/86400 * 24);
  
    % Surface temperature: Constant value TMEAN if TRANGE = 0. Otherwise, use a sine
    % wave with max (TMEAN + 1/2 TRANGE) at 2 pm and min (TMEAN - 1/2 TRANGE) at 2 am
    tsurf = tmean + 0.5 * trange * sin(2*pi/24 * (hour-8.0)); % 采用正弦函数来反映温度变化
    
    % Thermal conductivity and heat capacity
    soilvar = soil_thermal_properties (physcon, soilvar);
    
    % Soil temperatures
    soilvar = soil_temperature (physcon, soilvar, tsurf, dt);
    
    % Save hourly soil temperature profile for last day
    if (iday == nday)
      % Surface output
      % Vector format - to write to data file
      m = m + 1;
      hour_vec(m) = hour;
      z_vec(m) = 0;
      tsoi_vec(m) = tsurf - physcon.tfrz; % deg C
      
      % For MATLAB contour
      hour_out(itim) = hour;
      z_out(1) = 0;
      tsoi_out(1,itim) = tsurf - physcon.tfrz; % deg C
      
      % Soil layers for top 100 cm
      for i = 1:soilvar.nsoi
        if (soilvar.z(i) > -1.0)
          m = m + 1;
          hour_vec(m) = hour;
          z_vec(m) = soilvar.z(i) * 100; % cm
          tsoi_vec(m) = soilvar.tsoi(i) - physcon.tfrz; % deg C
          
          z_out(i+1) = soilvar.z(i) * 100; % cm
          tsoi_out(i+1,itim) = soilvar.tsoi(i) - physcon.tfrz; % deg C
        end
      end
    end
  end
end

% --- Write output file
A = [hour_vec; z_vec; tsoi_vec];
fileID = fopen('data.txt','w');
fprintf(fileID,'%12s %12s %12s\n','hour','z','tsoi');
fprintf(fileID,'%12.3f %12.3f %12.3f\n', A);
fclose(fileID);

% --- Make contour plot
contour(hour_out, z_out, tsoi_out,'ShowText','on')
title('Soil Temperature (^oC)')
xlabel('Time of day (hours)')
ylabel('Soil depth (cm)')
