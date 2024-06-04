% Supplemental program 7.1
clear,clc
import module_Radiation.*
import module_HydroTools.*
import module_Ipaper.*

sp_forcing

%% 30天之后的稳态
for j = 1:nday
  fprintf('day = %6.0f\n',j)
  
  for i = 1:ntim
    hour = i * (dt/86400 * 24); % Hour of day (0 to 24)
    
    % Air temperature (K): use a sine wave with max (tmean + 1/2 trange) at 1400
    % and min (tmean - 1/2 trange) at 0200
    forcvar.tref = tmean + 0.5 * trange * sin(2*pi/24 * (hour-8)) + physcon.tfrz; % Ta_ref
    
    % Vapor pressure (Pa) using constant relative humidity
    [esat, desat] = satvap (forcvar.tref-physcon.tfrz);
    forcvar.eref = (relhum / 100) * esat;
    
    % Derived quantities
    % forcvar.thref   ! Potential temperature at reference height (K)
    % forcvar.qref    ! Specific humidity at reference height (kg/kg)
    % forcvar.thvref  ! Virtual potential temperature at reference height (K)
    % forcvar.rhomol  ! Molar density at reference height (mol/m3)
    % forcvar.rhoair  ! Air density at reference height (kg/m3)
    % forcvar.mmair   ! Molecular mass of air at reference height (kg/mol)
    % forcvar.cpair   ! Specific heat of air at constant pressure, at reference height (J/mol/K)
    [forcvar.cpair, forcvar.thref, forcvar.thvref, forcvar.mmair, forcvar.rhomol] = ...
      cal_Cp(forcvar.tref, forcvar.eref, forcvar.pref, forcvar.zref);
    
    % Solar radiation (W/m2)
    % doy        ! Day of year (1 to 365)
    % lat        ! Latitude (radians)
    % decl       ! Declination (radians): Brock, T.D. (1981) Calculating solar radiation
    %            ! for ecological studies. Ecological Modelling 14:1-19
    % hour_angle ! Solar hour angle (radians)
    % coszen     ! Cosine of solar zenith angle
    % rv         ! Radius vector: Brock, T.D. 1981. Calculating solar radiation
    %            ! for ecological studies. Ecological Modelling 14:1-19
    % Rs_toa     ! Solar radiation on horizontal surface at top of atmosphere (W/m2)
    % tau_atm    ! Atmospheric transmission coefficient
    % oam        ! Optical air mass
    % soldir     ! Direct beam solar radiation on horizontal surface (W/m2)
    % soldif     ! Diffuse solar radiation on horizontal surface (W/m2)
    % solrad     ! Total solar radiation on horizontal surface (W/m2)
    
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
    
    % Calculate the soil temperatures and surface fluxes
    [fluxvar, soilvar, bucket] = surface_fluxes(physcon, forcvar, surfvar, soilvar, fluxvar, bucket, dt);
    
    % Rainfall to equal evaporative loss (kg H2O/m2/s)
    %     forcvar.rain = fluxvar.etflx * physcon.mmh2o;
    forcvar.rain = 0;
    
    % Bucket model hydrology
    [bucket] = bucket_hydrology (physcon, forcvar, fluxvar, bucket, dt);
    
    % Save data for graphics
    xhour(i) = hour;
    ytsrf(i) = fluxvar.tsrf - physcon.tfrz;
    ytref(i) = forcvar.tref - physcon.tfrz;
    yrnet(i) = fluxvar.rnet;
    yshflx(i) = fluxvar.shflx;
    ylhflx(i) = fluxvar.lhflx;
    ygsoi(i) = fluxvar.gsoi;
    yustar(i) = fluxvar.ustar;
    ygac(i) = fluxvar.gac * 100;
    ygcan(i) = surfvar.gcan * 100;
  end
end

% --- Write output files
A = [xhour; ytref; ytsrf; yrnet; yshflx; ylhflx; ygsoi; yustar; ygac; ygcan];
fileID = fopen('flux.txt','w');
fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n','hour','Ta','Ts','Rn','H','LE','G','ustar','gac','gcan');
fprintf(fileID,'%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n', A);
fclose(fileID);

B = [soilvar.z; soilvar.tsoi];
fileID = fopen('tsoi.txt','w');
fprintf(fileID,'%12s %12s\n','depth','tsoi');
fprintf(fileID,'%12.3f %12.3f\n', B);
fclose(fileID);

% --- Make graph
plot(xhour,yrnet,'g-',xhour,yshflx,'r-',xhour,ylhflx,'b-',xhour,ygsoi,'r--',xhour,ygac,'m-')
axis([0 24 -100 600])
set(gca,'xTick',0:3:24)
set(gca,'yTick',-100:100:800)
title('Diurnal cycle')
xlabel('Time of day (hours)')
ylabel('Flux (W m^{-2})')
legend('R_n','H','\lambdaE','G','g_{ac}*100','Location','northwest')
grid on;

% saveas(gcf,'Figure7_Surface_Energy_Fluxes.png')
