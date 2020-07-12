function [uE0, lb, info] = get_ue0(p, eb, lb0)
  %% created at 2007/07/27 by Bas Kooijman; modified 2010/05/02
  %% p: 1 or 3 -vector with parameters g, k_J/ k_M, v_H^b, see get_lb
  %% eb: optional scalar with scaled reserve density at birth (default eb = 1)
  %% lb0: optional scalar with scaled length at birth
  %%     default: lb is obtained from get_lb
  %% uE0: scalar with scaled reserve at t=0: U_E^0 g^2 k_M^3/ v^2
  %%      U_E^0 = M_E^0/ {J_EAm}
  %% lb: scalar with scaled length at birth l_b = L_b/ L_m
  %% info: scalar for failure (0) or success (1) of convergence

  if exist('eb', 'var') == 0
    eb = 1; % maximum value as juvenile
  end

  if exist('lb0', 'var') == 0
    if length(p) < 3
      fprintf('not enough input parameters, see get_lb \n');
      uE0 = []; lb = []; info = 0;
      return;
    end
    [lb, info] = get_lb(p, eb);
  else
    lb = lb0; info = 1;
  end

  %% unpack p
  g = p(1);  % energy investment ratio

  xb = g/ (eb + g);
  uE0 = (3 * g/ (3 * g * xb^(1/ 3)/ lb - beta0(0, xb)))^3;
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
