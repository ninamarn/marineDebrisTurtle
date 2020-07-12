function [pars, varargout] = pars_Caretta_caretta()
% ===== ★ ★ ★ ★ ★ =====
% pars_Caretta_caretta.m
% ©2020 Marko Jusup
% <mjusup[at]gmail.com>
% ===== ★ ★ ★ ★ ★ =====
% Function that specifies the DEB parameter values for loggerhead turtles.
% Also balances the macrochemical reactions that are essential to
% estimate respiration or metabolic water production. This balancing
% is currently implemented for illustrative purposes only.
% == energy-budget parameters ==
p_Am   = 747.271;                                                   % surface-area-specific assimilation rate (J/cm^2/d)
v          = 0.068074;                                                 % energy conductance (cm/d)
kap      = 0.72856;                                                   % allocation fraction to soma
p_M     = 11.2005;                                                   % volume-specific somatic maintenance rate (J/cm^3/d)
p_T      = 0.0;                                                           % surface-area-specific somatic maintenance rate (J/cm^2/d)
k_J       = 0.0011233;                                               % maturity maintenance rate coefficient (1/d)
E_G      = 7322;                                                        % volume-specific cost of structure (J/cm^3)
E_Hb    = 25350;                                                      % maturity at birth (J)
E_Hp    = 9.875e07;                                                 % maturity at puberty (J)
pars.ene=[p_Am, v, kap, p_M, p_T, k_J, E_G, E_Hb, E_Hp]';
% == chemical potentials ==
%               X       V       E       P
mu_O = [525, 500,  550,  480]' * 1000;                 % column vector of chemical potentials for organic compounds (J/Cmol)
% == chemical composition ==
%              X       V       E       P
n_O = [1.00, 1.00, 1.00, 1.00;                                % C/C, equals 1 by definition
            1.80, 1.80, 1.80, 1.80;                                % H/C, these values show that we consider dry-mass
            0.50, 0.50, 0.50, 0.50;                                % O/C
            0.15, 0.15, 0.15, 0.15];                               % N/C
%          C     H    O    N
n_M = [ 1     0     0     1;                                           % C/C, equals 0 or 1
             0     2     0     4;                                           % H/C
             2     1     2     1;                                           % O/C
             0     0     0     2];                                         % N/C
w_O = n_O' * [12, 1, 16, 14]';                                 % column-vector of molar dry masses for organic compounds (g/Cmol)
%            X      V      E      P
d_O = [0.3;  0.3;  0.3;  0.3];                                   % dry-to-wet mass ratio for organic compounds
M_V = d_O(2) / w_O(2);                                          % volume-specific amount of structure (Cmol/cm^3)
% assuming structural density of 1 g/cm^3
% == transformation efficiencies ==
kap_X=0.8;                                                              % fraction of ingested energy retained during assimilation
kap_G = M_V * mu_O(2) / E_G;                               % fraction of growth energy that gets built into structure
kap_P=0.1*kap_X;                                                   % energy in feces per unit of assimilated energy 
% note: kapP/kapX is the fraction of ingested energy lost to feces,
%          i.e., what journal articles call digestibility
kap_R=0.95;                                                            % fraction of reproduction energy turned into ova
pars.eff=[kap_X, kap_G, kap_P, kap_R]';
% == macrochemical reactions ==
% Assimilation
% yXE=1+cA1+cA4+yPE
% 1.80*yXE=1.80+2*cA2+4*cA4+1.80*yPE
% 0.50*yXE+2*cA3=0.50+2*cA1+cA2+cA4+0.50*yPE
% 0.15*yXE=0.15+2*cA4+0.15*yPE
% Growth
% 1=yVE+cG1+cG4
% 1.50=1.80*yVE+2*cG2+4*cG4
% 0.50+2*cG3=0.50*yVE+2*cG1+cG2+cG4
% 0.15=0.15*yVE+2*cG4
% Dissipation
% 1=cD1+cD4
% 1.50=2*cD2+4*cD4
% 0.50+2*cD3=2*cD1+cD2+cD4
% 0.15=2*cD4
% Macrochemical reactions with four elements have three degrees of freedom.
% These are the yields yXE, yPE, and yVE. To fix them, three efficiencies must be defined.
% The three efficiencies are kapX, kapP, and kapG, where
% yXE=kapX^(-1)*muE/muX
% yPE=kapP*muE/muP
% yVE=kapG*muE/muV
y_XE=kap_X^(-1)*mu_O(3)/mu_O(1);
y_PE=kap_P*mu_O(3)/mu_O(4);
y_VE=kap_G*mu_O(3)/mu_O(2);
c_A=n_M\(n_O(:,1)*y_XE-n_O(:,3)-n_O(:,4)*y_PE);
c_G=n_M\(n_O(:,3)-n_O(:,2)*y_VE);
c_D=n_M\n_O(:,3);
bal.yie=[y_XE, y_PE, y_VE]';
bal.A=c_A;
bal.G=c_G;
bal.D=c_D;
varargout{1}=bal;
% == auxiliary parameters ==
E_Hh      = 21080;                                                   % maturity at hatching if different from EHb
E_Hj       = NaN;                                                       % maturity at metamorphosis if needed
T_A        = 7200;                                                     % Arrhenius temperature
del_M     = 0.38997;                                                % Shape factor
pars.aux = [E_Hh, E_Hj, T_A, del_M]';
% == useful compound parameters ==
k_M = p_M / E_G;                                                     % somatic maintenance rate coefficient (1/d)
k = k_J / k_M;                                                           % maintenance ratio
E_m = p_Am / v;                                                      % reserve capacity (J/cm^3)
g = E_G / kap / E_m;                                                % energy investment ratio
L_m = v / k_M / g;                                                    % maximum structural length (cm)
L_T = p_T / p_M;                                                      % heating length (cm)
l_T = L_T / L_m;                                                       % scaled heating length
pars.com = [k_M, k, E_m, g, L_m, L_T, l_T]';
% == other useful conversions ==
rho_E = 1000 * mu_O(3) * d_O(3) / w_O(3);           % mass-reserve coupler (J/kg)
E2u = k_M / p_Am / L_m^2;                                    % energy 2 unergy, i.e., scaled energy (1/J)
rho_U = rho_E * E2u;                                               % mass-unergy coupler (1/kg)
rho_F = 1000 * mu_O(1)*d_O(1) / w_O(1) * E2u;   % food-unergy coupler (1/kg)
pars.con = [rho_E, E2u, rho_U, rho_F]';
% == scaled threshold maturities ==
u_Hh = E_Hh * E2u;                                                 % scaled maturity at hatching
u_Hb = E_Hb * E2u;                                                 % scaled maturity at birth
u_Hj = E_Hj * E2u;                                                   % scaled maturity at metamorphosis
u_Hp = E_Hp * E2u;                                                 % scaled maturity at puberty
pars.thr = [u_Hh, u_Hb, u_Hj, u_Hp]';
return
