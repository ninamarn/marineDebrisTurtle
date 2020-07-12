% ===== ★ ★ ★ ★ ★ =====
% pop_Caretta_caretta_wrapper.m
% ©2020 Marko Jusup
% <mjusup[at]gmail.com>
% ===== ★ ★ ★ ★ ★ =====
% This script is a wrapper for running 'pop_Caretta_caretta.m', which sets up a projection matrix
% to calculate the population growth rate of loggerhead turtles given the following environmental
% conditions:
% A) default food abundance;
% B) default temperature;
% 1) percentage of gastrointestinal tract occupied by plastic debris;
% 2) residence time of debris relative to food.
% One simulation is run for each set of conditions 1) and 2), which then modulate the default
% value of A). Condition B) is fixed throughout. If temperature is changed, it is also necessary
% to take into account the temperature-dependent sex determination in loggerhead turtles.
% The results are stored in the matrix 'lammat'. This matrix was used to generate fig. 4 in the study
% (but with many more points; settings here are for illustration to keep the run time reasonable).
more off                                                                   % let Octave print freely to screen
clear('-x', 'lammat', 'apmat', 'elmat')                      % cleanup vars from previous run if needed
global Food Temp
Food0 = 0.809;                                                       % set default food abundance
Temp = 21.8;                                                          % set default temperature
% == settings ==
pctV1=linspace(0.01,0.03,3);                                 % percentage of gastrointestinal tract occupied by debris
pctV2=linspace(0.03,0.25,20);
pctV3=linspace(0.25,0.50,20);
pctV=[pctV1(1:end), pctV2(2:end), pctV3(2:end)];
RelResTime=linspace(1,10,10);                             % residence time of debris relative to food
% resolution is increased by increasing the
% third parameter of the linspace function
% == space allocation ==
Nrows=length(pctV);
Ncols=length(RelResTime);
if ~exist('lam', 'var')
  lammat = zeros(Nrows, Ncols);                            % for population growth rate
  apmat = zeros(Nrows, Ncols);                              % for age of sexual maturation
  elmat = cell(Nrows, Ncols);                                  % elasticity matrices (not used in the study)
end;
% == simulations ==
for r = 1 : Nrows
  for c = 1 : Ncols
    RYX=RelResTime(c) * pctV(r) / ( 1-pctV(r) );     % rel. residence time multiplied by debris-to-food ratio
    Food = Food0 / (1+Food0 * RYX);                       % implement eq. [3] from the study
    pop_Caretta_caretta;                                          % build projection matrix based on current conditions
    apmat(r, c) = ap / 365.25;                                  % record age to sexual maturation in years
    [lammat(r, c), elmat{r, c}] = get_lam(Amat);  % calculate population growth and elasticities
    % == displaying status ==
    disp(['row=' num2str(r) ', col=' num2str(c) ', lam=' num2str(lammat(r, c))])
    disp(' ')
    pause(1)
  end;
end;