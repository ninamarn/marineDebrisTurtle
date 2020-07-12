function out = get_F(u_H)
    % ===== ★ ★ ★ ★ ★ =====
    % get_F.m
    % ©2020 Marko Jusup
    % <mjusup[at]gmail.com>
    % ===== ★ ★ ★ ★ ★ =====
    % calculates reproductive output
    % ...
    % when u_H >= u_Hp, energy invested into reproduction
    % between time moments t and t+1 is u_H(t+1)-u_H(t)
    % one egg costs u_E0, but conversion efficiency must not be forgotten
    global pars
    kap_R  = pars.eff(4);
    u_Hp   = pars.thr(4);
    global u_E0
    out = (u_H(2:end)>=u_Hp)*kap_R.*diff(u_H)/u_E0;
end