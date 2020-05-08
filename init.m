function simu = init(fs,ff)
%---------------------------------------------------------------
% Define simulation and individual features
%
% simu: structure with individual features (parameters, env...)
%           (simu(i).par = parameters of individual i )
%
% called by:main.m
% calls: set_par.m
%
% created: 2013/03/12 - Laure Pecquerie
% modified: 2018/10/30 - Nina Marn
%---------------------------------------------------------------

% setup
 t0 = 1; % initial time (hatching /start of feeding) % 1st July (midseason)
 % originally August 1st --> important when having seasonal T !!! 
 tm = 67 * 365; % final time (duration of simulation)
 ts = 60; % hatching time (same as spawning ?)
  
 % color for plots
% red = [1 0 0] ; %  30   % -50%
% orange = [1 0.4 0]; % 28  %  
% yello = [1 1 0]; % 26  %   -30%
% green = [0 0.6 0]; % 24  % -20%
% blue = [0 0 1]; % 22  % - 10%
% dblue= [0 0 0.4]; % 20  %  0
% purple = [0.6 0 0.6]; % 18 % +20%
% black = [0 0 0 ]; % 16  % +50%

env = 1; % 1 = constant, 2 = fluctuating
 

style1 = '-';
style2 = '--';
% style = [repmat({style1},[length(fs)-1,1]); {style2}];
% sty = style{ff}
sty = style1; 

% color = [0 0.749019622802734 0.749019622802734]; % light blue
% color = [0.749019622802734 0.749019622802734 0 ]; % yellow-green
% color = {'b';'k'};
% colormap = parula; 
cmap=parula;
cmap(:,1)=linspace(min(cmap(:,1)),max(cmap(:,1)),64)';
% cmap(:,2)=linspace(min(cmap(:,2)),max(cmap(:,2)),64)'*0.5;
% colormap(cmap)
jump = 9 - 3 ;  % number of simulations-1
ind = [1: jump : 64, 64]; % index for plotting --> form 64 rows of colormap, choose steps
col = cmap(ind(ff),:);
% col = cmap(ff,:);
markers = repmat({'*' ; 'd' ; 'v' ; 'o'; '^'},[length(fs),1]); 
mark = markers{ff}; % marker default
lgdTxt = num2str(fs);

 % initial forcing variable
T = 21.8 + 273; 
% fs = [1, 0.81, 0.75, 0.656, 0.613]; % 1- max, 0.81 - in nature; 0.751 = 3% debris, KX=KY , 
% % 0.656 = 3% debris, KX = 3KY; 0.613 = 25% debris KX=KY; 
 f = fs(ff);
 
 % initiate state variable values
 V_0 = 0000.1; % cm^3, initial structural volume
 E_H0 = 0; % J, initial cum.energy invested into maturation
 E_R0 = 0; % J, reproduction buffer

%%% set parameters
  simu.par = set_par; % load parameters
%%%%%
 
 
% simulation 
simu.t0 = t0; 
simu.tm = tm; 
simu.ts = ts; 
simu.env = env; 
simu.Tinit = T; 
simu.finit = f;  % f = x/ (1+x) = X/K / (1 + X/K)


% % individual
p_Am = simu.par(5);
v = simu.par(6);
kap = simu.par(7);
p_M = simu.par(9);
kJ = simu.par(11);
E_G = simu.par(12);
E_Hb = simu.par(13);
VHb = (E_Hb/p_Am)/ (1 - kap);
g = (E_G/kap) / (p_Am/v);
kM = p_M /E_G;

p_UE0 = [VHb; g; kJ; kM; v]; % pars for initial scaled reserve
[UE0 Lb info] = initial_scaled_reserve (f, p_UE0);
E_0 = UE0 * p_Am;

simu.EVHR_init = [E_0 ; V_0 ; E_H0 ; E_R0]; % initial vector for individual
simu.col = col; 
simu.sty = sty; 
simu.mark = mark;
simu.lgdTxt = lgdTxt;
