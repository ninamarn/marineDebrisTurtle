function [d, ps] = getELHR_starv(t, ELHR, cPar, par, TC_env, tFs)

% 2019-09-19 Nina Marn; calculates individual fluxes while including starvation rules
% called by new_simu.m
% made for the marine debris ingestion paper; modeling work by N Marn and M Jusup

if length(tFs)>1
    F = spline1(t, tFs);
    if F >0.809
        fprintf('Oh NO, f = %2.2.f \n', F)
    elseif F <0.710
        fprintf('Hungry times, f = %2.2f \n', F)
    end
else
    F = tFs;
end

%% unpack state vars
E  = ELHR(1); % J, reserve energy
L  = ELHR(2); % cm, structural length
E_H  = ELHR(3); % J , cumulated energy inversted into maturity
E_R  = ELHR(4); % J, reproduction buffer


% unpack par, data, auxData
vars_pull(cPar); vars_pull(par);
pT_Am  = p_Am * TC_env; vT = v*TC_env; kT_M = k_M * TC_env; pT_M = p_M* TC_env;
kT_J = k_J *TC_env;

pA = pT_Am*F*L^2; % assimilation flux
% fprintf(1, 'e = %2.2f.  \n ', F, (E./L.^3)/E_m);
pC = (E_m *(vT/L + kT_M) * F * g/(F+g))*L^3; % mobilization flux (eq.2.20 without L_T)
pM = pT_M*L^3; % somatic maintenance
pG = kap*pC - pM; % growth
pJ = kT_J * E_H; % maturity maintenance
pR = (1-kap)*pC - pJ; % to maturation / reproduction

 % Differential equations
if  (pC*kap >= pM) || pC*kap > pM * 0.995 %(abs(pC*kap - pM) < 1e4) % we assume (pM + pJ) < pC % if enough mobilized to satisfy the maintenance, carry on as usual
     ps = [pA; pC; pM; pG; pJ; pR]; 
      dE = pA - pC; % dE/dt
      if E_H < E_Hp
          dH = pR; % dEH/dt
          dR = 0; % dER/dt
      else
          dH = 0;
          dR = pR;
      end

else % kap * pC < pM % energy needed to satisfy maintenance
     del_pC = -pG; % or del_pC = pC - (pM + pJ); % how much is missing for maintenance
     pG = 0; % if pM > kap*pC , then prevent shrinking, and remember the difference
     if E_R > 0 % if there is something in the repro buffer
         pR = (1-kap)*pC - pJ; % take the missing energy from repro buffer
         pC_eff = pC + del_pC;
         dH = 0;
         dR = pR - del_pC; 
         dE = pA - pC; % dE/dt
         fprintf(1, 'hungry, I took %2.4f kJ from ER, can still make %4.2f eggs! \n', del_pC/1e3, E_R/1e3/210);
     else
         %          keyboard
         if kap*pC >= pM * 0.998;
            fprintf(1,'This is idiotic !!!! \n')
             return
         end
         pC_eff = pM + pJ;
         fprintf(1, 'starving, no E_R, while f = %2.2f, and e = %2.2f.  \n ', F, (E./L.^3)/E_m);
         fprintf(1, 'Btw, pC is currently %2.2f MJ, razlika kap*pC i pM = %2.2f MJ  \n ', pC/1e6, abs(kap*pC/1e6 - pM/1e6) );
         dE = pA - pC_eff; pR = 0; % uzima samo osnovno iz rezerve + crisis, no eggs for you
         %       dE = pA - pC; pR = pC - pM - pJ; % povlaèi sve što može
         % first satisfy somatic maintenance, then whatever is left goes to maturity maintenance
         dH = 0; % this could result in rejuvenation if pR<0 (ie pC-pM < pJ)
         dR = pR;
     end
     %        ps = [pA; pC; pM; pG; pJ; pR; pC_eff; pM+pJ ] / 1e06
     
end

  dL=1/(3*L^2)*pG/ E_G;  %dL/dt
  %dV = (kap * pC - pS) / E_G;% dV/dt
 
  d = [dE; dL; dH; dR];
  
    
%   disp(dE)
%   disp(dL)
%   disp(dH)
%   disp(dR)

end
