
function d = flux(t, EVHR, simu)
  %---------------------------------------------------------------
  % Define differential equations of the state variables
  %  
  % t: n-vector with time points
  % EVHR: 4-vector with state variables
  %         E , J, reserve energy
  %         V , cm^3, structural volume
  %         E_H , J , cumulated energy inversted into maturity (E_H in Kooijman 2010)
  %         E_R , J, reproduction buffer (E_R in Kooijman 2010)
  %         
  % d: 4-vector with d/dt E, V, E_H, E_R
  %
  % called by : indiv.m, 
  % calls : food.m, temp.m, 
  %
  % 2013/03/15 - Laure Pecquerie
  %--------------------------------------------------------------
  % Environmental conditions
  if simu.env == 1
      f = simu.finit ; % keep food constant
      T = simu.Tinit;  % temp constant
  else
      f = scaled_f_resp(t);
      T = simu.Tinit;  % temp constant
  end
   
  %% unpack state vars
  E  = EVHR(1); % J, reserve energy
  V  = EVHR(2); % cm^3, structural volume
  E_H  = EVHR(3); % J , cumulated energy inversted into maturity 
  E_R  = EVHR(4); % J, reproduction buffer
  
  % read parameter values
  par = simu.par;
  T1 = par(1);    % K, Reference temperature ; 
  TA = par(2);    % K, Arrhenius temperature ;
  kap_X = par(4); % -, digestion efficiency of food to reserve
  p_Am = par(5);  % J/cm^2/d, maximum surface-specific assimilation rate
  v = par(6);      % cm/d, energy conductance
  kap = par(7);    % -, allocation fraction to soma = growth + somatic maintenance
  p_M = par(9);    % J/d.cm^3, [p_M], vol-specific somatic maintenance
  p_T = par(10);   % J/d.cm^2, {p_T}, surface-specific som maintenance
  k_J = par(11);   % 1/d, maturity maint rate coefficient
  E_G = par(12);   % J/cm^3, [E_G], spec cost for structure
  E_Hb = par(13);  % J, E_H^b, maturity at birth
  E_Hp = par(14);  % J, E_H^p, maturity at puberty
  K = par(22);     % same unit as food, half saturation constant
  
    
  % simplest temperature correction function, 
  % see Kooijman 2010 for more detailed formulation (e.g. p. 21)
  c_T = exp(TA/ T1 - TA/ T) ;
  p_AmT = c_T * p_Am ;
  v_T = c_T * v; % 
  p_MT = c_T * p_M;
  p_TT = c_T * p_T; % 
  k_JT = c_T * k_J; 
  p_XmT = p_AmT / kap_X;
  
  % Fluxes
  if E_H < E_Hb
      pX = 0;% embryo stage
  else
      pX = f * p_XmT * V^(2/3);
  end
  
  pA = kap_X * pX;
  pM = p_MT * V;
  pT = p_TT * V^(2/3);
  pS = pM + pT;
  pC = (E/V) * (E_G * v_T * V^(2/3) + pS ) / (kap * E/V + E_G ); %eq. 2.12 p.37 Kooijman 2010
  pJ = k_JT * E_H;
  
  % Differential equations
  dE = pA - pC; % dE/dt
  dV = (kap * pC - pS) / E_G;% dV/dt
  if E_H < E_Hp
      dH = (1 - kap) * pC - pJ; % dEH/dt
      dR = 0; % dER/dt
  else
      dH = 0;
      dR = (1 - kap) * pC - pJ;
  end
  
  d = [dE; dV; dH; dR]; 