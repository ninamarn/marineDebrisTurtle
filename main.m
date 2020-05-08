%---------------------------------------------------------------
% Compute model predictions for one individual,
%   with standard DEB equations (Kooijman 2010)
% 
% Parameter estimates for Caretta caretta
%
% from birth (first feeding), to 10 weeks of age (SCL) / or 4 years
% (weight)
% no reproduction
%(shrinking of structural volume is possible)
% 
% simu : structure with individual features and model predictions
%
% calls :   init.m
%           indiv.m (calls flux.m)
%           get_obs.m, get_plots2 
%
% 2018/10/30, 2019/10/29 - Nina Marn, modified from scripts by Laure Pecquerie
% 
%---------------------------------------------------------------  

  clear all
  close all
  
  %% 1 - initialize time, parameters, etc
%   fs = fliplr( linspace(0.81, 0.613+0.002,10)); % 1- max, 0.81 - in nature; 0.751 = 3% debris, KX=KY , 
    % 0.656 = 3% debris, KX = 3KY; 0.613 = 25% debris KX=KY; 
fs = fliplr( linspace(0.81, 0.64 + 0.002,10)); % 1- max, 0.81 - in nature; 0.64 at digestive system occupancy of 25%
    
  for ff = 1 : length(fs)
       simu = init(fs, ff) ;
  
  %% 2 - calculate state variables
  simu.tEVHR = indiv(simu);
  
  %% 3 - calculate observable quantities
  simu.obs = get_obs(simu);
  
  %% 4 - make plots
  hold on
    get_plots2(simu)
  
  end
  
  return
  XX = 11;
  YY = 10;
  fNo = 0; 
  
  %%% add legend
  figure(1+fNo)
  %    subplot(3,1, 1)
%   legend(simu.lgdTxt)
  set(gcf,'PaperPosition',[0 0 XX YY]);
  print('-r200','-dtiff','Fig_DEB_a.tif');
  
  
  figure(2+fNo)
  set(gcf,'PaperPosition',[0 0 XX YY]);
  print('-r300','-dtiff','Fig_DEB_b.tif');
  
  figure(3+fNo)
  set(gcf,'PaperPosition',[0 0 XX+2 YY+2]);
  print('-r300','-dtiff','Fig_DEB_c.tif');
  
  
  