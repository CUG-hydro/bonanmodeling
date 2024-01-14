classdef module_Radiation
  %% module_Radiation
  %
  %% Example:
  % import module_Radiation.*
  
  properties (Constant)
    solcon = 1364; % W/m2
  end
  
  properties
  end
  
  methods (Static)
    function hello()
      fprintf("hello world\n");
    end
    
    function [Rs_toa, Rs, Rs_dir, Rs_dif, coszen] = cal_Rs_toa(lat, doy, hour)
      % Solar radiation at top of the atmosphere
      solcon = 1364; % W/m2
      
      decl = 23.45 * sin((284+doy)/365*2*pi) * pi/180;
      hour_angle = 15 * (hour-12) * pi/180;
      coszen = max(cos(lat)*cos(decl)*cos(hour_angle) + sin(lat)*sin(decl), 0);
      rv = 1 / sqrt(1 + 0.033*cos(doy/365*2*pi));
      Rs_toa = solcon / rv^2 * coszen;
      
      % Clear sky atmospheric attenuation: Gates, D.M. (1980) Biophysical Ecology, page 110, 115
      tau_atm = 0.5;
      oam = 1 / max(coszen, 0.04);
      Rs_dir = Rs_toa * tau_atm^oam;                       % Clear sky direct beam
      Rs_dif = Rs_toa * (0.271 - 0.294 * tau_atm^oam);     % Clear sky diffuse
      Rs = Rs_dif + Rs_dir;                                % Clear sky total
    end
    
    function [cpair, thref, thvref, mmair, rhomol] = cal_Cp(tref, eref, pref, zref)
      % thref   ! Potential temperature at reference height (K)
      % qref    ! Specific humidity at reference height (kg/kg)
      % thvref  ! Virtual potential temperature at reference height (K)
      % rhomol  ! Molar density at reference height (mol/m3)
      % rhoair  ! Air density at reference height (kg/m3)
      % mmair   ! Molecular mass of air at reference height (kg/mol)
      % cpair   ! Specific heat of air at constant pressure, at reference height (J/mol/K)
      thref =  tref + 0.0098 * zref; % T_pot
      qref = physcon.mmh2o / physcon.mmdry * eref / ...
        (pref - (1 - physcon.mmh2o/physcon.mmdry) * eref);
      
      thvref = thref * (1 + 0.61 * qref); % Tv
      
      rhomol = pref / (physcon.rgas * tref);
      rhoair = rhomol * physcon.mmdry * ...
        (1 - (1 - physcon.mmh2o / physcon.mmdry) * eref / pref);
      mmair = rhoair / rhomol;
      
      cpair = physcon.cpd * (1 + (physcon.cpw/physcon.cpd - 1) * qref); % J/kg/K
      cpair = cpair * mmair;                                            % J/mol/K
    end
  end
end

