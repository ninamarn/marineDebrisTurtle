function out = get_W(a, u_E, l, u_H)
    % ===== ★ ★ ★ ★ ★ =====
    % get_W.m
    % ©2020 Marko Jusup
    % <mjusup[at]gmail.com>
    % ===== ★ ★ ★ ★ ★ =====
    % calculates body mass
    % ...
    % uses the average state of the reproductive buffer in the adult stage
    global pars
    L_m    = pars.com(5);
    rho_U  = pars.con(3);
    u_Hp   = pars.thr(4);
    a_p = interp1(u_H, a, u_Hp, 'pchip');                  % age at sexual maturation
    dela = 0 + (u_H>=u_Hp) .* min(a-a_p, .5);       % time passed since sexual maturation
    u_R = u_H - interp1(a, u_H, a-dela, 'pchip');      % unergy in reproductive buffer
    out = (1.0e-3 * (l*L_m).^3 + (u_E + u_R) / rho_U); % body mass...
    % consists of three components:
    % 1) the mass of structure
    % 2) the mass of reserve
    % 3) the mass of reproductive buffer
end