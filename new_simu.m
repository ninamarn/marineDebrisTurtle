% 2019-09-19 Nina Marn; simulate some fluctuating debris 

clear all 
close all
%% Environment
envTemp = C2K(21.8);  % Hawkes et al. 2011
f = 0.81;  % Marn et al., 2017

%% define %V to simulate: V_Y/V_X , values based on Frick et al 2009
Frick_mean = 0.03; % real value: 3.4%  --> percent of TOTAL content! (so, V_Y / (V_Y + V_X) !!! correct in calculations)
Frick_max = 0.25; % real value: 25.7%
sim_max = 0.50; % max in simulation

% make a simple t-Vs vector, time is in days
time = [0 : 30 :  64*365 ]; % step of 1 month, up to cca 65 years


% r = pearsrnd(mu,sigma,skew,kurt,m,n) returns an m-by-n matrix of random numbers drawn 
%from the distribution in the Pearson system with mean mu, standard deviation sigma, skewness skew, and kurtosis kurt. 
% The parameters mu, sigma, skew, and kurt must be scalars. Normal -> kurt = 3
% mu = Frick_mean;
% mu = 0.1; % future scenario
% sigma = 0.04;
% skew = 2.5  + 0.3; 
% kurt = 70 + 25 + 10;
m = 1; 
n = length(time)*2;
% %  r = pearsrnd(mu,sigma,skew,kurt,m,n);
% % save('new_simu', 'r')
%  r_new = pearsrnd(mu,sigma,skew,kurt,m,n);
% save('new_simu', 'r_new', '-append')
%  load new_simu

% set distribution - lognormal
mu = 0.6; sigma = 1; 
x = 0:0.1:30; y = pdf('lognormal', x, mu, sigma );  %--> to visualize
% plot(x,y, 'k-', 'LineWidth', 2)

r = lognrnd(mu, sigma, m,n); % --> pull numbers
r(r>30) = []; % remove all of those where % plastics in GI would be >30%

Vs = r(1:length(time)); 
%  Vs = 0; % without plastic
 load new_simu_plastics_19-OCt-2019.mat % <-- this should override the Vs with old ones
 
% Verify the values in r are within the specified range.
%  V_range = [min(Vs) max(Vs)]

Y2X = Vs/100 ./ (1-Vs/100);


%% load parameters, setup the organism etc
load results_Caretta_caretta.mat
cPar = parscomp_st(par);
vars_pull(cPar); vars_pull(par)

% temp corrections
 birthTemp = C2K(29); % temp experienced until birth; see Marn et al. 2017 (DEB logg)
 pubTemp = envTemp; % temp experienced until puberty is lower (see Marn et al. 2017 (DEB logg)) but for simplicity 21.8 assumed 

 TC_env = tempcorr(envTemp, T_ref, T_A);
 TC_ab = tempcorr(birthTemp, T_ref, T_A);
 TC_tp = tempcorr(pubTemp, T_ref, T_A);
 
% feeding & assimilation
FT_m = F_m * TC_env;       % 3, l/d.cm^2, {F_m} max spec searching rate
pT_Am = p_Am * TC_env;       % 5, J/cm^2/d, maximum surface-specific assimilation rate
K_X = pT_Am /(kap_X*FT_m); 

%% run simulation
X = f*K_X / (1-f);
thetas = Y2X; %  omjer theta_Y i theta_X = omjer V_Y i V_X ; % staro: thetas = (Vs).^(2/3); % row vector 
Ys = thetas * X; % row vector

K = K_X*(1+ Ys/K_X );
fs = X ./ (X+K); % <-- resulting fs , from that debris load

 % initial state
  pars_tp = [g k l_T v_Hb v_Hp];  % life cycle
  [tau_p, tau_b, l_p, l_b, info] = get_tp(pars_tp, f);
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  ELHR_b = [f*E_m  L_b E_Hb 0]; % initial state
 
 % puberty and ultimate
 L_p = L_m * l_p;                  % cm, structural length at puberty at f
 L_i = L_m * f;                  % cm, ultimate structural length at f
 
 % physical length and age
 Lw_b = L_b/del_SCL; Lw_p = L_p/del_SCL; Lw_i = L_i/del_SCL
 aT_b = (tau_b/k_M + t_0)/ TC_ab;    % d, age at birth at f and T
 tT_p = (tau_p/k_M)/ TC_tp + aT_b;   % d, age at puberty at f and T
 age_at_pub081 = tT_p/365
 
 %%% 
% %   tc = ts; % set current time to initial time
%   ELHR_tc = ELHR_b; % set current values of state variables to initial values
%  tELHR = zeros(0,5); % initialize outputs
%   i = 0; % initialize day index

 % solve for state variables
 if Vs == 0
     tFs = f;
 else
     tFs = [time' , fs'];
 end
 
 ELHR_tc = ELHR_b ; 
 tELHR = zeros(0,5); % initialize outputs
 i = 0; iR = 0; iS = 0; % initialize year index, and repro index, and simulation index
 ts = 0; tc = ts; 
   while tc < time(end)
      i = i + 2; % nesting every two years
      tnext = ts + 365*i; % integrate two nesting dates (one nesting each year)
      if tnext > time(end)
          tnext  = time(end);
      end
      
      tsimu = [tc; tnext]';
      [t, ELHR, ps] = ode45(@(tsimu, ELHR)getELHR_starv(tsimu, ELHR, cPar, par, TC_env, tFs), tsimu, ELHR_tc); % integrate
%   [t, EVHR] = ode45(flux(t, EVHR, simu), t, EVHR_tc) % integrate
      tELHR = [tELHR; [t,ELHR]]; % append to output
      E_Hc = ELHR(end,3);
      if E_Hc >= E_Hp % E_Hp
          iR = iR+1; ER = ELHR(end,4);
          trepro(iR) = tsimu(end); eggs(iR) = floor(ER/210/1e3); 
          fprintf(1,'Reproducing, go my %2.2f eggs! \n', eggs(iR));
           ELHR(end,4) = ER - eggs(iR)*210*1e3; % reproduce
      end
      
      ELHR_tc = ELHR(end, :)'; tc = tsimu(end);
  end     
 

% [time , ELHR , ps] = ode45(@(tsimu, ELHR)getELHR_starv(tsimu, ELHR, cPar, par, TC_env, tFs ), tsimu, ELHR_b);
E =  tELHR(:,2); L =  tELHR(:,3); E_H =  tELHR(:,4); R =  tELHR(:,5);
tsimed = tELHR(:,1);

pT_Am  = cPar.p_Am * TC_env; vT = v*TC_env; kT_M = k_M * TC_env; pT_M = p_M* TC_env; kT_J = k_J *TC_env; % correct for temp
% read life history
 Lw = L/del_SCL; 
ind_Lwi = find(Lw >= (0.995*Lw_i), 1); % at this point, 99.5% of ultimate length reached
aT_Li = tsimed(ind_Lwi); aT_Li/365;
ind_pub = find(R > 0, 1); 
Lwp = Lw(ind_pub); aTp = tsimed(ind_pub)/365 

% total body mass = structure + (reserve + reprobuffer)
%W_V = d_V * L.^3; % dry weight of structure
%W_E_and_R = w_E / mu_E * (E + R);  % dry W of E and E_R
%W = W_V + W_E_and_R; % total dry weight
%Ww = 1/d_V * W / 1000;% kg assume the same water content in structure and reserves

rho_E = w_E/mu_E/d_V
Ww = L.^3 + rho_E *(E + R) /1000; % kg assume the same water content in structure and reserves

return
if Vs == 0
    normal_Lw = Lw; normal_Ww = Ww; normal_t = tsimed; normal_eggs = eggs; normal_trepro = trepro;
    norm_E = E; norm_E_H = E_H; norm_R = R; norm_L = L;
    save('new_simu_normal.mat', 'normal_Lw', 'normal_Ww', 'normal_t', 'normal_eggs', 'normal_trepro', ...
        'norm_E','norm_E_H', 'norm_R', 'norm_L')
else 
    eval(['save(''new_simu_plastics_', date ,...
        '.mat'', ''Lw'', ''Ww'', ''tsimed'', ''eggs'', ''trepro'', ''E'', ''R'', ''E_H'', ''L'', ''tsimed'', ''Vs'', ''tFs'')'])
end

%% plotting
load new_simu_normal
normal_BF = 100 * normal_Ww*1000 ./ normal_Lw.^3; 
BF = 100 * Ww*1000 ./ Lw.^3; 
[counts, binCntrs] = hist(Vs,50);
hold on
n = length(Vs);
binWidth = binCntrs(3) - binCntrs(2); 
prob = counts / (n * binWidth);


figure
subplot(2,2,1)
hold on
H = bar(binCntrs,prob,'hist');
set(H,'facecolor',[0.5 0.5 0.5]);
plot(x,y, 'k-', 'LineWidth', 2)
legend('normalized %debris in dig.cont.', 'pdf of the assumed %debris distribution')
title('simulated %debris in dig.cont.')

subplot(2,2,2)
plotyy(time/365, Vs, time/365, fs)
title('% debris in dig.cont and resulting f')

subplot(2,2,3)
hold on
plot(normal_t/365, normal_BF, 'b')
plot(tsimed/365, BF, 'r')
title('Body condition factor over time')

subplot(2,2,4)
hold on
plot(normal_trepro/365, normal_eggs, 'bo')
plot(trepro/365, eggs, 'ro')
title('Reproduction over time')

return  % extra figures

figure
subplot(2,2,1)
hold on
plot(normal_t/365, normal_Lw, 'b')
plot(tsimed/365, Lw, 'r')
title('Length in time')

subplot(2,2,2)
hold on
plot(tsimed/365, Ww, 'r')
plot(normal_t/365, normal_Ww, 'b')
title('Mass in time')

% plote state variables
% clear all
% % load new_simu_plastics_17-Oct-2019
load new_simu_normal
figure
subplot(2,2,1)
hold on
title('Structural length over time')
plot(normal_t, norm_L ,'b')
plot(tsimed, L ,'r')

subplot(2,2,2)
hold on
title('Maturity over time')
plot(normal_t, norm_E_H ,'b')
plot(tsimed, E_H ,'r')

subplot(2,2,3)
hold on
title('Reserve over time')
plot(normal_t, norm_E ,'b')
plot(tsimed, E,'r')

subplot(2,2,4)
hold on
title('Repro buffer over time')
plot(normal_t, norm_R ,'b')
plot(tsimed, R ,'r')

figure
hold on
title('e=[E]/[E_m] over time')
plot(normal_t/365, (norm_E./norm_L.^3)/E_m)
plot(tsimed/365, (E./L.^3)/E_m)


% saving figures:
return 
  print('-r300','-dtiff','fig_newSimu_Oct2019_main_detail.tif');
  print('-r300','-dtiff','fig_newSimu_Oct2019_L_Ww.tif');
  