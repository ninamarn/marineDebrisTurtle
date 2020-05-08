function get_plots2(simu)

  par = simu.par; % parameters
  E_Hp = par(14);
  K = par(22); 
  Em = par(5) / par(6); % J/cm^3  [Em] = pAm / v

  % time , environmental forcing
  t = simu.tEVHR(:,1) / 365;
  
  if simu.env == 1 
          f = repmat(simu.finit, 1, length( t ) );
          T = repmat(simu.Tinit, 1, length( t) ) - 273; % in DegC
  else
          f = scaled_f_resp(t);
          T = repmat(simu.Tinit, 1, length( t) ) - 273; % in DegC
  end
  
    
  % state variables
  E = simu.tEVHR(:,2);
  V = simu.tEVHR(:,3);
  E_H = simu.tEVHR(:,4);
  E_R = simu.tEVHR(:,5);
  
  sc_res_dens = E ./ ( V * Em); % -, scaled reserve density e
  
  % spawning date indices
  i_sp = find(and((E_R == 0),(E_H>=E_Hp)));
  i_sp = i_sp - 1 ; % look at the preceding line with E_R value before spawning
  
  %markup 
  color = simu.col; 
%   sty = simu.sty;
%   marker = simu.mark ;
 marker = 'none' ;
  dots = [1 : 365*8 : t(end)*365]; 
  dotsF = [1: 2 : length(i_sp)];
%   t = t(dots); T = T(dots); f = f(dots);
  
  % observable quantities
  Lw = simu.obs(:,2); 
  Ww = simu.obs(:,3)/1000; 
  Ew = simu.obs(:,4); 
  F  = simu.obs(:,5); 
      
  % graph
%   figure
%   set(gcf,'PaperPositionMode','manual');
%   set(gcf,'PaperUnits','centimeters');
  %left bottom width height
% %   set(gcf,'PaperPosition',[0 0 8.5 25]); 
  
%   subplot(2, 3, 1)
%   hold on
%   plot( t, T, color, 'Linewidth',1.5)
%   plot( t(dots), T(dots), [color marker] )
%   title('Temperature')
%   xlabel('Time (yr)')
%   ylabel('T (degC)')
%   axis([0 t(end) 0 30])

%   subplot(2,3, 2)
%   hold on
%   plot( t(dots), f(dots), [color marker '-'], 'Linewidth',1.5)
%   title('Scaled functional response')
%   xlabel('Time (yr)')	
%   ylabel( 'f ')
%   axis([0 t(end) 0 1.2])
  
fNo = 0; 
  % Lw
  figure(1+fNo)
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','centimeters');
%   subplot(3,1, 1)
  hold on
%   keyboard
  plot(t, Lw, 'Color', color, 'Linewidth',1.5)
%   plot(t(dots), Lw(dots), 'Color', color, 'Marker', marker, 'Linewidth',1.5)
  title('A')
  xlabel('Time (yr)')
  ylabel(' Carapace length (cm)')
  axis([0 t(end)+2 0 max(Lw)+10])
  
  % Ww_
figure(2+fNo)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
  %   subplot(3,1, 2)
  hold on
  plot( t,  Ww, 'Color', color, 'Linewidth',1.5)
%   plot( t(dots),  Ww(dots), 'Color', color, 'Marker', marker, 'Linewidth',1.5)
  title('B')
  xlabel('Time (yr)')
  ylabel('Body mass (kg)')
  axis([0 t(end)+2 0 max(Ww)+10])
  
  % F
  figure(3+fNo)
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','centimeters');
%   subplot(3, 1, 3)
  hold on
  plot(t(i_sp), F(i_sp), 'Color', color, 'Marker', '.', 'Linewidth',1.5, 'LineStyle', 'none', 'MarkerSize', 20)
%   plot(t(i_sp(dotsF)), F(i_sp(dotsF)), 'Color', color, 'Marker', marker, 'Linewidth',1.5)
  title('C')
  xlabel('Time (yr)')
  ylabel('Fecundity per nesting (#)')
  axis([0 t(end)+2 0 max(F)+20])
  % Create colorbar
% colorbar('peer',axes1,'LineWidth',1);

% F
  figure(4+fNo)
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','centimeters');
  hold on
  plot(Lw, Ww, 'Color', color, 'Linewidth',1.5)
  title('D')
  xlabel('Carapace length (cm)')
  ylabel('Body mass (kg)')
  axis([0 max(Lw)+10 , 0 max(Ww)+10])
  % Create colorbar
% colorbar('peer',axes1,'LineWidth',1);


     
%   print -djpeg90 fig_DEB_standard4.jpg 

