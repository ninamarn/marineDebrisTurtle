function [lb, info] = get_lb(p, eb, lb0)
  %% [lb, info] = get_lb(p, eb, lb0)
  %% created 2007/08/15 by Bas Kooijman; modified 2011/07/27
  %% p: 3-vector with parameters: g, k_J/ k_M, v_H^b (see below)
  %% eb: optional scalar with scaled reserve density at birth (default eb = 1)
  %% lb0: optional scalar with initial estimate for scaled length at birth
  %% lb: scalar with scaled length at birth
  %% info: indicator equals 1 if successful
  %% solves y(x_b) = y_b  for lb with explicit solution for y(x)
  %%   y(x) = x e_H/(1-kap) = x g u_H/ l^3
  %%   and y_b = x_b g u_H^b/ ((1-kap)l_b^3)
  %%   d/dx y = r(x) - y s(x);
  %%    with solution y(x) = v(x) \int r(x)/ v(x) dx
  %%    and v(x) = exp(- \int s(x) dx)
  %% replacement of Euler integration by ode23: get_lb1
  %%    get_lb is much faster than get_lb1
  %% shooting method: get_lb2

  %% unpack p
  g = p(1);  % energy investment ratio: 
  k = p(2);  % k_J/ k_M, ratio of maturity and somatic maintenance rate coeff
  vHb = p(3); % v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_B^b = M_H^b/ {J_EAm}
  info = 1;

  if exist('lb0','var') == 0
    lb0 = [];
  end
   
  if k == 1
    lb = vHb^(1/ 3); % exact solution for k = 1
    info = 1;
    return;
  end
  if isempty(lb0)
    lb = vHb^(1/ 3); % exact solution for k = 1     
  else
    lb = lb0;
  end
  if exist('eb','var') == 0
    eb = 1;
  elseif isempty(eb)
    eb = 1;
  end
     
  n = 1000 + round(1000 * max(0, k - 1)); xb = g/ (g + eb); xb3 = xb ^ (1/3);
  x = linspace(1e-6, xb, n); dx = xb/ n;  x3 = x .^ (1/3);
  
  b = beta0(x, xb)/ (3 * g); t0 = xb * g * vHb;
  i = 0; norm = 1; % make sure that we start Newton Raphson procedure
  ni = 100; % max number of iterations
  
  while i < ni  && norm > 1e-8
    l = x3 ./ (xb3/ lb - b);
    s = (k - x) ./ (1 - x) .* l/ g ./ x;
    v = exp( - dx * cumsum(s)); vb = v(n);
    r = (g + l); rv = r ./ v;
    t = t0/ lb^3/ vb - dx * sum(rv);
    dl = xb3/ lb^2 * l .^ 2 ./ x3; dlnl = dl ./ l;
    dv = v .* exp( - dx * cumsum(s .* dlnl));
    dvb = dv(n); dlnv = dv ./ v; dlnvb = dlnv(n);
    dr = dl; dlnr = dr ./ r;
    dt = - t0/ lb^3/ vb * (3/ lb + dlnvb) - dx * sum((dlnr - dlnv) .* rv);
    %% [i lb t dt] % print progress
    lb = lb - t/ dt; % Newton Raphson step
    norm = t^2; i = i + 1;
  end
    
  if i == ni || lb < 0 || lb > 1 || isnan(norm) % no convergence
    % try to recover with a shooting method
    if isempty(lb0)
      [lb, info] = get_lb2(p, eb);
    elseif lb0 < 1 && lb0 > 0
      [lb, info] = get_lb2(p, eb, lb0);
    else
      [lb, info] = get_lb2(p, eb);
    end
  end
  
  if info == 0
    fprintf('warning get_lb: no convergence of l_b \n');
  end
end

function [lb info] = get_lb2(p, eb, lb0)
  %% [lb info] = get_lb2(p, eb, lb0)
  %% created 2007/07/26 by Bas Kooijman; modified 2013/08/21
  %% p: 4-vector with parameters: g, k_J/ k_M, v_H^b (see below)
  %% eb: optional scalar with scaled reserve density at birth (default eb = 1)
  %% lb0: optional scalar with initial estimate for scaled length at birth
  %% lb: scalar with scaled length at birth
  %% info: indicator equals 1 if successful
  %% solves y(x_b) = y_b  with shooting method
  %%   y(x) = x e_H = x g u_H/ l^3 and y_b = x_b g u_H^b/ l_b^3
  %% calls dget_lb2 and fnget_lb2

  global xb xb3 kap g k vHb
  
  %% unpack p
  g = p(1);   % energy investment ratio
  k = p(2);   % k_J/ k_M, ratio of maturity and somatic maintenance rate coeff
  vHb = p(3); % v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_B^b = M_H^b/ {J_EAm}

  info = 1;

  if exist('lb0','var') == 0 || k == 1
    lb = vHb ^ (1/3); % exact solution for k = 1
    if k == 1
      return;
    end
  elseif isempty(lb0)
    lb = vHb^(1/ 3); % exact solution for k = 1     
  else
    lb = lb0;
  end
  if exist('eb','var') == 0
    eb = 1;
  elseif isempty(eb)
    eb = 1;
  end

  xb = g/ (eb + g); xb3 = xb^(1/ 3);
  
  options = optimset('Display','off');
  [lb, flag, info] = fzero('fnget_lb2',lb,options);
  if lb < 0 || lb > 1
      info = 0;
  end
end

function f = fnget_lb2(lb)
  %% f = y(x_b) - y_b = 0; x = g/ (e + g); x_b = g/ (e_b + g)
  %% y(x) = x e_H = x g u_H/ l^3 and y_b = x_b g u_H^b/ l_b^3

  global xb g vHb Lb

  Lb = lb;
  y = lsode('dget_lb2', 0, [1e-10, xb]);
  f = y(2) - xb * g * vHb/ lb^3;
end

function dy = dget_lb2(y,x)
  %% y = x e_H/ (1 - kap); x = g/ (g + e)
  %% (x,y): (0,0) -> (xb, xb eHb/ (1 - kap)) 

  global Lb xb xb3 g k

  x3 = x^(1/ 3);
  l = x3/ (xb3/ Lb - beta0(x, xb)/ 3/ g);
  dy = l + g - y * (k - x)/ (1 - x) * l/ g/ x;
end

function f = beta0 (x0,x1)
  %% f = beta0 (x0,x1)
  %% created 2000/08/16 by Bas Kooijman, modified 2011/04/10
  %% special incomplete beta function:
  %%   B_x1(4/3,0) - B_x0(4/3,0) = \int_x0^x1 t^(4/3-1) (1-t)^(-1) dt

  if x0 < 0 | x0 >= 1 | x1 < 0 | x1 >= 1
    fprintf('Warning from beta0: argument values outside (0,1) \n');
    f = [];
    return;
  end

  n0 = length(x0); n1 = length(x1);
  if n0 ~= n1 && n0 ~= 1 && n1 ~= 1
    fprintf('Warning from beta0: argument sizes don not match \n');
    f = [];
    return;
  end
  
  x03 = x0 .^ (1/ 3); x13 = x1 .^ (1/ 3); a3 = sqrt(3);
  
  f1 = - 3 * x13 + a3 * atan((1 + 2 * x13)/ a3) - log(x13 - 1) + ...
      log(1 + x13 + x13 .^ 2)/ 2;

  f0 = - 3 * x03 + a3 * atan((1 + 2 * x03)/ a3) - log(x03 - 1) + ...
      log(1 + x03 + x03 .^ 2)/ 2;
  
  f = f1 - f0;
end
