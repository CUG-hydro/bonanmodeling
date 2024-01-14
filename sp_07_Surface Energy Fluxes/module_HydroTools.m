classdef module_HydroTools
  %UNTITLED 此处提供此类的摘要
  %   此处提供详细说明
  properties (Constant)
  end
  
  methods (Static)    
    function [phi_c] = phi_c_monin_obukhov (x)
      % --- Evaluate the Monin-Obukhov phi function for scalars at x
      if (x < 0)
        phi_c = (1 - 16 * x)^(-0.5);
      else
        phi_c = 1 + 5 * x;
      end
    end
    
    function [phi_m] = phi_m_monin_obukhov (x)
      % --- Evaluate the Monin-Obukhov phi function for momentum at x
      if (x < 0)
        phi_m = (1 - 16 * x)^(-0.25);
      else
        phi_m = 1 + 5 * x;
      end
    end
    
    function [psi_c] = psi_c_monin_obukhov (x)
      % --- Evaluate the Monin-Obukhov psi function for scalars at x
      if (x < 0)
        y = (1 - 16 * x)^0.25;
        psi_c = 2 * log((1 + y^2)/2);
      else
        psi_c = -5 * x;
      end
    end
    
    function [psi_m] = psi_m_monin_obukhov (x)
      % --- Evaluate the Monin-Obukhov psi function for momentum at x
      if (x < 0)
        y = (1 - 16 * x)^0.25;
        psi_m = 2 * log((1 + y)/2) + log((1 + y^2)/2) - 2 * atan(y) + pi / 2;
      else
        psi_m = -5 * x;
      end
    end
    
    function [psi_hat_c] = psi_c_rsl (z, h, L, c1, c2)
      % --- Evaluate the roughness sublayer (RSL) function psi_hat for scalars
      % at z. Note that z has already been adjusted for the displacement height
      % (i.e., using z - d).
      % ------------------------------------------------------
      % Input
      %   z            ! Vertical height - displacement height (m)
      %   h            ! Canopy height - displacement height (m)
      %   L            ! Obukhov length (m)
      %   c1           ! Parameter for RSL function phi_hat (dimensionless)
      %   c2           ! Parameter for RSL function phi_hat (dimensionless)
      %
      % Output
      %   psi_hat_c    ! RSL psi_hat function for scalars (dimensionless)
      % ------------------------------------------------------
      % The function to integrate depends on unstable (f1) or stable (f2)
      f1 = @(x) (1-16*x/L).^(-0.5) .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
      f2 = @(x) (1+5*x/L)          .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
      
      % Numerically integrate the function from z to infinity
      if (L < 0)
        psi_hat_c = integral (f1, z, inf);
      else
        psi_hat_c = integral (f2, z, inf);
      end
    end
    
    
    function [psi_hat_m] = psi_m_rsl (z, h, L, c1, c2)
      % --- Evaluate the roughness sublayer (RSL) function psi_hat for momentum
      % at z. Note that z has already been adjusted for the displacement height
      % (i.e., using z - d).
      
      % ------------------------------------------------------
      % Input
      %   z            ! Vertical height - displacement height (m)
      %   h            ! Canopy height - displacement height (m)
      %   L            ! Obukhov length (m)
      %   c1           ! Parameter for RSL function phi_hat (dimensionless)
      %   c2           ! Parameter for RSL function phi_hat (dimensionless)
      %
      % Output
      %   psi_hat_m    ! RSL psi_hat function for momentum (dimensionless)
      % ------------------------------------------------------
      
      % The function to integrate depends on unstable (f1) or stable (f2)
      f1 = @(x) (1-16*x/L).^(-0.25) .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
      f2 = @(x) (1+5*x/L)           .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
      
      % Numerically integrate the function from z to infinity
      if (L < 0)
        psi_hat_m = integral (f1, z, inf);
      else
        psi_hat_m = integral (f2, z, inf);
      end
    end
    
    
    function [fluxvar, fx] = most (physcon, forcvar, surfvar, fluxvar, x)
      
      % Use Monin-Obukhov similarity theory to obtain the Obukhov length (obu).
      % This is the function to solve for the Obukhov length. For the current
      % estimate of the Obukhov length (x), calculate u*, T*, and q* and then
      % the new length (obu). The function value is the change in Obukhov length:
      % fx = x - obu.
      
      % -------------------------------------------------------------------------
      % Input
      %   x                  ! Current estimate for Obukhov length (m)
      %   physcon.vkc        ! von Karman constant
      %   physcon.grav       ! Gravitational acceleration (m/s2)
      %   physcon.mmh2o      ! Molecular mass of water (kg/mol)
      %   forcvar.zref       ! Reference height (m)
      %   forcvar.uref       ! Wind speed at reference height (m/s)
      %   forcvar.thref      ! Potential temperature at reference height (K)
      %   forcvar.thvref     ! Virtual potential temperature at reference height (K)
      %   forcvar.eref       ! Vapor pressure at reference height (Pa)
      %   forcvar.pref       ! Atmospheric pressure (Pa)
      %   forcvar.mmair      ! Molecular mass of air at reference height (kg/mol)
      %   fluxvar.tsrf       ! Surface temperature (K)
      %   fluxvar.esrf       ! Surface vapor pressure (Pa)
      %   fluxvar.z0m        ! Roughness length for momentum (m)
      %   fluxvar.z0c        ! Roughness length for scalars (m)
      %   fluxvar.disp       ! Displacement height (m)
      % Output
      %   fluxvar.ustar      ! Friction velocity (m/s)
      %   fluxvar.tstar      ! Temperature scale (K)
      %   fluxvar.qstar      ! Water vapor scale (mol/mol)
      %   fluxvar.obu        ! Obukhov length (m)
      %   fx                 ! Change in Obukhov length (x - obu)
      % -------------------------------------------------------------------------
      
      % --- Prevent near-zero values of the Obukhov length
      if (abs(x) <= 0.1)
        x = 0.1;
      end
      
      % --- Calculate z-d at the reference height, because this is used many times
      z_minus_d = forcvar.zref - fluxvar.disp;
      
      % --- Evaluate psi for momentum at the reference height (zref-disp) and surface (z0m)
      [psi_m_zref] = module_HydroTools.psi_m_monin_obukhov (z_minus_d / x);
      [psi_m_z0m] = module_HydroTools.psi_m_monin_obukhov (fluxvar.z0m / x);
      psim = -psi_m_zref + psi_m_z0m;
      
      % --- Evaluate psi for scalars at the reference height (zref-disp) and surface (z0c)
      [psi_c_zref] = module_HydroTools.psi_c_monin_obukhov (z_minus_d / x);
      [psi_c_z0c] = module_HydroTools.psi_c_monin_obukhov (fluxvar.z0c / x);
      psic = -psi_c_zref + psi_c_z0c;
      
      % --- Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)
      
      zlog_m = log(z_minus_d / fluxvar.z0m);
      zlog_c = log(z_minus_d / fluxvar.z0c);
      
      fluxvar.ustar = forcvar.uref * physcon.vkc / (zlog_m + psim);
      fluxvar.tstar = (forcvar.thref - fluxvar.tsrf) * physcon.vkc / (zlog_c + psic);
      fluxvar.qstar = (forcvar.eref - fluxvar.esrf) / forcvar.pref * physcon.vkc / (zlog_c + psic);
      tvstar = fluxvar.tstar + 0.61 * forcvar.thref * fluxvar.qstar * (physcon.mmh2o / forcvar.mmair);
      
      % --- Calculate Obukhov length (m)
      fluxvar.obu = fluxvar.ustar^2 * forcvar.thvref / (physcon.vkc * physcon.grav * tvstar);
      
      % --- Calculate change in obu
      fx = x - fluxvar.obu;
    end
    
    function [esat, desat] = satvap (tc)
      % Compute saturation vapor pressure and change in saturation vapor pressure
      % with respect to temperature. Polynomial approximations are from:
      % Flatau et al. (1992) Polynomial fits to saturation vapor pressure.
      % Journal of Applied Meteorology 31:1507-1513. Input temperature is Celsius.
      
      % --- For water vapor (temperature range is 0C to 100C)
      a0 =  6.11213476;        b0 =  0.444017302;
      a1 =  0.444007856;       b1 =  0.286064092e-01;
      a2 =  0.143064234e-01;   b2 =  0.794683137e-03;
      a3 =  0.264461437e-03;   b3 =  0.121211669e-04;
      a4 =  0.305903558e-05;   b4 =  0.103354611e-06;
      a5 =  0.196237241e-07;   b5 =  0.404125005e-09;
      a6 =  0.892344772e-10;   b6 = -0.788037859e-12;
      a7 = -0.373208410e-12;   b7 = -0.114596802e-13;
      a8 =  0.209339997e-15;   b8 =  0.381294516e-16;
      
      % --- For ice (temperature range is -75C to 0C)
      c0 =  6.11123516;        d0 =  0.503277922;
      c1 =  0.503109514;       d1 =  0.377289173e-01;
      c2 =  0.188369801e-01;   d2 =  0.126801703e-02;
      c3 =  0.420547422e-03;   d3 =  0.249468427e-04;
      c4 =  0.614396778e-05;   d4 =  0.313703411e-06;
      c5 =  0.602780717e-07;   d5 =  0.257180651e-08;
      c6 =  0.387940929e-09;   d6 =  0.133268878e-10;
      c7 =  0.149436277e-11;   d7 =  0.394116744e-13;
      c8 =  0.262655803e-14;   d8 =  0.498070196e-16;
      
      % --- Limit temperature to -75C to 100C
      tc = min(tc, 100);
      tc = max(tc, -75);
      
      % --- Saturation vapor pressure (esat, mb) and derivative (desat, mb)
      if (tc >= 0)
        esat  = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 ...
          + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8)))))));
        desat = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 ...
          + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8)))))));
      else
        esat  = c0 + tc*(c1 + tc*(c2 + tc*(c3 + tc*(c4 ...
          + tc*(c5 + tc*(c6 + tc*(c7 + tc*c8)))))));
        desat = d0 + tc*(d1 + tc*(d2 + tc*(d3 + tc*(d4 ...
          + tc*(d5 + tc*(d6 + tc*(d7 + tc*d8)))))));
      end
      
      % --- Convert from mb to Pa
      esat  = esat  * 100;
      desat = desat * 100;
    end
    
    % function obj = untitled(inputArg1,inputArg2)
    %   obj.Property1 = inputArg1 + inputArg2;
    % end
    % function outputArg = method1(obj,inputArg)
    %   outputArg = obj.Property1 + inputArg;
    % end
  end
end
