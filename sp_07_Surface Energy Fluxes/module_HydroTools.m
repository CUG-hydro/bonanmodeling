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
    
    
    function [fluxvar, fx] = most (x, varargin)
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
      if length(varargin) == 1 && iscell(varargin{1}); varargin = varargin{1}; end
      [physcon, forcvar, ~, fluxvar] = flatten(varargin);
      % [physcon, forcvar, surfvar, fluxvar] = flatten(varargin);
      
      % --- Prevent near-zero values of the Obukhov length
      if (abs(x) <= 0.1); x = 0.1; end
      
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
    
  end
end
