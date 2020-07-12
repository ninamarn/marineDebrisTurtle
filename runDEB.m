function varargout = runDEB(y0, a, varargin)
    % ===== ★ ★ ★ ★ ★ =====
    % runDEB.m
    % ©2020 Marko Jusup
    % <mjusup[at]gmail.com>
    % ===== ★ ★ ★ ★ ★ =====
    % Function accompanied by a small set of auxiliary functions (see below)
    % that together run the DEB model. Parameter values for any organism
    % can be supplemented (see file 'pars_Caretta_caretta.m' for an example
    % how to do this). The DEB model is cast in a scaled form using the following
    % definitions:
    % r_*  = p_* / ( {p_Am} × (l×L_m)^2 )                 (energy flows)
    % u_* = E_* × k_M / ( {p_Am} × L_m^2 )            (energy)
    % l = L / L_m                                                         (structural length)
    % This form of the model is convenient for a couple of reasons:
    % 1) state variables are dimensionless;
    % 2) initial energy reserve in terms of u_E is finite;
    % 3) utilisation flow r_C is non-zero when l=0.
    % Note that state-dependent or time-varying food and temperature can be
    % defined in the files 'f_ya.' and 'T_ya.m'.
    global pars bal                                                    % secure access to parameter values and mass balance
    % == energy-budget parameters ==
    global kap                                                          % make parameter values accessible to subroutines
    p_Am   = pars.ene(1);                                         % surface-area-specific maximum assimilation rate (J//cm^2/d)
    v          = pars.ene(2);                                         % energy conductance (cm/d)
    kap      = pars.ene(3);                                         % allocation fraction to soma
    p_M      = pars.ene(4);                                        % volume-specific somatic maintenance rate (J/cm^3/d)
    p_T       = pars.ene(5);                                        % surface-area-specific somatic maintenance rate (J/cm^2/d)
    k_J        = pars.ene(6);                                        % maturity maintenance rate coefficient (1/d)
    E_G      = pars.ene(7);                                        % volume-specific cost of structure (J/cm^3)
    E_Hb    = pars.ene(8);                                        % maturity at birth (J)
    E_Hp    = pars.ene(9);                                         % maturity at puberty (J)
    % == efficiencies ==
    global kap_X kap_G kap_R                                  % make parameter values accessible to subroutines
    kap_X  = pars.eff(1);                                           % assimilation efficiency
    kap_G  = pars.eff(2);                                           % growth efficiency
    kap_P  = pars.eff(3);                                           % egestion efficiency
    kap_R  = pars.eff(4);                                           % reproduction efficiency
    % == auxiliary parameters ==
    global T_A                                                           % make parameter values accessible to subroutines
    E_Hh    = pars.aux(1);                                         % maturity at hatching (J), EHh=EHb as a first approximation
    E_Hj     = pars.aux(2);                                         % maturity at metamorphosis (J)
    T_A      = pars.aux(3);                                         % Arrhenius temperature (K)
    del_M   = pars.aux(4);                                         % shape coefficient
    % == compound parameters ==
    global k_M k g l_T                                               % make parameter values accessible to subroutines
    k_M      = pars.com(1);
    k          = pars.com(2);
    E_m     = pars.com(3);
    g          = pars.com(4);
    L_m      = pars.com(5);
    L_T       = pars.com(6);
    l_T        = pars.com(7);
    % == scaled threshold maturities ==
    global u_Hb u_Hj u_Hp                                       % make parameter values accessible to subroutines
    u_Hh    = pars.thr(1);
    u_Hb    = pars.thr(2);
    u_Hj     = pars.thr(3);
    u_Hp    = pars.thr(4);
    % == prepare ==
    global u_E0 l_b l_j                                               % make quantities accessible to subroutines
    % initial conditions are in y0
    % state vector is y=[u_EA, u_EC, l, u_H]
    if isempty(y0)                                                     % if initial conditions are unspecified we need to specify them now
        f0 = f_ya();                                                     % mother's food abundance
        l_b = get_lb([g, k, u_Hb/(1-kap)], f0);            % scaled length at birth
        u_E0 = get_ue0(g, f0, l_b);                             % scaled initial energy reserve
        u_Eb = f0 * l_b^3 / g;                                      % scaled energy reserve at birth
        y0 = [u_E0, u_E0-u_Eb, l_b, u_Hb];
    elseif length(y0) ~= 4
        error('I would appreciate the intial conditions for 4 state variables!')
    end;
    if isempty(a)
        error('How about supplying some information on the simulation time span?')
    elseif length(a) == 1
        a = linspace(0, a, a+1)';
    else
        a = linspace(a(1), a(end), a(end)-a(1)+1)';    
    end;
    opts = odeset('RelTol',1e-6,'AbsTol',1e-12);
    % == calculate length at metamorphosis ==
    l_j = l_b;                                                               % do NOT delete; must be here
    if ~isnan(u_Hj)
        l_j = inf;                                                           % do NOT delete; must be here
        [~, y] = ode45(@dydu_H, [u_Hb, u_Hj]', [0 y0(1:3)], opts);
        l_j = y(end,4); % length-at-metamorphosis
    end;
    % == model run ==
    [a, y] = ode45(@dyda, a, y0, opts);
    varargout{1} = y;
    varargout{2} = a;
    clear('-global', 'kap')
    clear('-global', 'kap_X', 'kap_G', 'kap_R')
    clear('-global', 'T_A')
    clear('-global', 'k_M', 'k', 'g', 'l_T')
    clear('-global', 'u_Hb', 'u_Hj', 'u_Hp')
end
function dy = dydu_H(u_H, y)
    % function to calculate length-at-metamorphosis
    % scaled maturity is the independent variable here
    % variables are analogous to the function below
    global kap
    global k_M k u_Hp
    global l_b l_j
    dy = zeros(4, 1);
    a = y(1); u_E = y(2)-y(3); l = y(4);
    [r_A, r_C, r_G] = flows(a, [y(2:4); u_H]);
    p_R = (1 - kap) * l^2 * r_C - k * min(u_H, u_Hp);
    T = T_ya();
    k_MT = k_M * C_T(T);
    dy(1) = 1 / (k_MT * p_R);
    dy(2) = l^2 * r_A / p_R;
    dy(3) = l^2 * r_C / p_R;
    dy(4) = r_G / 3 / kap / p_R;
end
function dy = dyda(a, y)
    % the model; this is where important things happen
    % age is the independent variable here
    % p_* = {p_Am} * (l*L_m)^2 * r_*
    global kap
    global k_M
    dy = zeros(4, 1);
    u_E = y(1) - y(2); l = y(3); u_H = y(4);
    [r_A, r_C, r_G, r_R] = flows(a, y);
    T = T_ya();
    k_MT = k_M * C_T(T);                                          % somatic maintenance rate coefficient needs to be adjusted to ambient temperature
    dy(1) = k_MT * l^2 * r_A;                                    % model eq. 1; assimilated energy
    dy(2) = k_MT * l^2 * r_C;                                    % model eq. 2; utilized energy
    dy(3) = k_MT * r_G / 3 / kap;                               % model eq. 3; growth
    dy(4) = k_MT * l^2 * r_R;                                    % model eq. 4; maturation
end
function varargout = flows(a, y)
    global kap
    global kap_X kap_G kap_R
    global k g l_T
    global u_Hb u_Hj u_Hp
    global l_b l_j
    u_E = y(1) - y(2); l = y(3); u_H = y(4);
    M1 = max(1, min(l/l_b, l_j/l_b));
    f_a = 0 + (u_H >= u_Hb) * f_ya(l, a);                % just in case one wants to play with time dependent food availability
    r_A = M1 * f_a;                                                    % scaled assimilation flow
    r_C = u_E * (M1 * g + l_T + l) / (l^3 + u_E);      % scaled utilization (mobilization) flow
    r_S = kap * (l + l_T);                                           % scaled somatic maintenance flow
    kap2 = max(kap, r_S/r_C);                                 % this is inserted to handle starvation reasonably
    r_G = kap2 * r_C - r_S;                                        % scaled growth flow
    r_J = k * max(0, min(u_H, u_Hp)/l^2);               % scaled maturity maintenance flow
    r_R = (1 - kap2) * r_C - r_J;                                 % scaled maturation (reproduction) flow
    % == gather outputs ==
    varargout{1} = r_A;
    varargout{2} = r_C;
    varargout{3} = r_G;
    varargout{4} = r_R;
end
function out = C_T(T)
    % temperature correction function
    global T_A
    T0 = 293.15;                                                       % reference temperature
    out = exp(T_A / T0 - T_A ./ T);
end