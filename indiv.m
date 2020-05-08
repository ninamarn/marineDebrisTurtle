function tEVHR = indiv(simu)
  %---------------------------------------------------------------
  % Compute model predictions (numerical integration)
  %  
  % simu: 1-structure with individual features (parameters, env...),;
  %             see init.m
  %
  % tEVHR : (nt,4) matrix with time and state variables
  %         t, d, time
  %         E , J, reserve energy
  %         V , cm^3, structural volume
  %         E_H , J , cumulated energy inversted into maturity
  %         E_R , J, reproduction buffer
  %         
  % called by : main.m
  % calls : fluxes.m
  %
  % 2013/03/15 - Laure Pecquerie
  %--------------------------------------------------------------
  

  par = simu.par;
  tc = simu.t0; % set current time to initial time
  EVHR_tc = simu.EVHR_init; % set current values of state variables to initial values
  
  tEVHR = zeros(0,5); % initialize outputs
  i = 0; % initialize year index
  while tc < simu.tm
      i = i + 2; % nesting every second year
      tnext = simu.ts + 365 * i; % integrate two spawning dates (one spawning each year)
      if tnext > simu.tm
          tnext = simu.tm;
      end
%    keyboard   
      t = [tc:tnext]'; 
      [t, EVHR] = ode45(@(t, EVHR)flux(t, EVHR, simu), t, EVHR_tc); % integrate
      tEVHR = [tEVHR; [t, EVHR]]; % append to output
      E_Hc = EVHR(end,3);
      if E_Hc >= par(14) %; E_Hp
          EVHR(end,4) = 0;
      end

      EVHR_tc = EVHR(end,:)'; tc = t(end);
      
  end   