clear all;close all
%%%  paper_Marine debris
% % A script for going from %volume to %AUs, and then to f
% Nina Marn, 2018/11/02; modified: 2019/07/26, 2019/09/26

%% Environment
T = C2K(21.8);  % Hawkes et al. 2011
f_0 = 0.81;  % Marn et al., 2017

%% define %V to simulate: V_Y/V_X , values based on Frick et al 2009
Frick_mean = 0.03; % real value: 3.4%  --> percent of TOTAL content! (so, V_Y / (V_Y + V_X) !!! correct in calculations)
Frick_max = 0.25; % real value: 25.7%

% debris-to-food ratio (in the ENVIRONMENT)
sim_max = 0.5; % max in simulation -- > 50% more plastics than food

Y2X = linspace(0,sim_max, 2000); 

R_K = [1; 1.25; 2; 3; 5; 10];% How much longer is residence time of debris relative to that of food
R_K2 = linspace(1,10, 1000); % for image surface
surface = 1; 
% R_K = K_X / K_Y, and a K is proportional to k_i (release time), and inv. proportional to T_i (res time) if b_i are equal

%  --> %V / (1- %V) = th_Y/th_X = K_X/K_Y * Y/X = R_K * Y2X
% --> %V = R_K*Y/X / (1+ R_K*Y/X) = Vs


%% run simulation

f_eff =  zeros(length(Y2X), length(R_K));
Vs =  zeros(length(Y2X), length(R_K));
for kk = 1 : length(R_K)
    Vs(:, kk ) = R_K(kk)*Y2X ./ (1+ R_K(kk)*Y2X) ; % each column is different residence time; rows are plastic conc
    f_eff(:, kk ) = f_0 ./ (1 + f_0 * R_K(kk) * Y2X);  
end

% Vs larger than 100% are impossible
iImposs = find(Vs>1);
Vs(iImposs) = NaN; f_eff(iImposs) = NaN; 

ind3 = find(0.025<=Vs & Vs<=0.035); % locate where %V = 3%
ind25 = find(0.24<=Vs & Vs<=0.26); % locate where %V = 25%

multiY2X = repmat(Y2X,1,length(R_K))'; % make a column of env debris-to-food ratios

if surface
   Vs2 =  zeros(length(Y2X), length(R_K2));
 for kk = 1 : length(R_K2)
    Vs2(:, kk ) = R_K2(kk)*Y2X ./ (1+ R_K2(kk)*Y2X)  ; % each column is different residence time; rows are plastic conc
 end
  [ind3r, ind3c]  = find(0.0298<=Vs2 & Vs2<=0.0302); % locate where %V = 3%
  [ind25r, ind25c] = find(0.2498<=Vs2 & Vs2<=0.2502); % locate where %V = 25%
  [ind50r, ind50c] = find(0.4998<=Vs2 & Vs2<=0.5002); % locate where %V = 50%
end



%% plot
LiWi = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]; % line Width -- darker line = slower release time
LiWi = LiWi *1.5; 

figure(1) % will make figure in manuscript
% subplot(2,1,1)
% hold on
% hold on
% for kk = 1: length(R_K)
% plot(Y2X, Vs(:,kk)*100, 'k', 'LineWidth', LiWi(kk) )
% end
% plot([0 sim_max], [3 3], 'k--')
% plot([0 sim_max], [25 25], 'k--')
% plot(multiY2X(ind3), Vs(ind3)*100, 'rx'  ) % za kontrolu
% plot(multiY2X(ind25), Vs(ind25)*100, 'bx'  ) % za kontrolu
% ylabel('% debris in digestive contents')

% subplot(2,1,2)
hold on
for kk = 1: length(R_K)
plot(Y2X, f_eff(:,kk), 'k', 'LineWidth', LiWi(kk) )
end
plot([min(Y2X) max(Y2X)], [mean(f_eff(ind3)) mean(f_eff(ind3))], 'k--'  ) % 
plot([min(Y2X) max(Y2X)],  [mean(f_eff(ind25)) mean(f_eff(ind25))], 'k--'  ) % 
xlabel('Environmental debris-to-food ratio')
ylabel('Perceived food abundance')

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
% print('-r300','-dtiff','debris2f_Oct2019_01.tif');
 


if surface
    figure
    oh=surf(R_K2, Y2X, Vs2);
    view(2)
    oh.EdgeColor='none';
    cmap = colormap(1); cmap2 = flipud(cmap); % get current colormap and modify
    colormap(1, cmap2) % change in figure to match other figures
    xlabel('Residence time of debris relative to food')
    ylabel('Environmental debris-to-food ratio')
    
%   print('-r300','-dtiff','debris2f_Oct2019_colors.tif');
 %  print('-r300','-dtiff','debris2f_Oct2019_colors.tif');
%     
    figure
    hold on
     for ii = 1: length(ind3r)
         row =ind3r(ii); clmn =  ind3c(ii);
         plot(R_K2(clmn), Y2X(row), 'k.', 'LineWidth', 1)
     end
     for ii = 1: length(ind25r)
         row =ind25r(ii); clmn =  ind25c(ii);
         plot(R_K2(clmn), Y2X(row), 'b.', 'LineWidth', 1)
     end
     for ii = 1: length(ind50r)
         row =ind50r(ii); clmn =  ind50c(ii);
         plot(R_K2(clmn), Y2X(row), 'b.', 'LineWidth', 1)
     end
     axis([1 10 0 sim_max])
     xlabel('Residence time of debris relative to food')
    ylabel('Environmental debris-to-food ratio')
    
%   print('-r300','-dtiff','debris2f_Oct2019_03.tif');
  
end


  
