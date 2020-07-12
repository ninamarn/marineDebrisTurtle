function out = get_S(a, l)
    % ===== ★ ★ ★ ★ ★ =====
    % get_S.m
    % ©2020 Marko Jusup
    % <mjusup[at]gmail.com>
    % ===== ★ ★ ★ ★ ★ =====
    % calculates survival from literature-specified hazard rates
    % ecologically, hazard is determined by size not age
    % returns survival up to all ages specified in 'a'
    % but for that purpose needs l=l(a) as produced by the DEB model
    global pars
    del_M = pars.aux(4);
    L_m = pars.com(5);
    l_i = del_M*[4 45 72 92]/L_m;
    h_i = -log([0.786 0.875 0.810 0.810 0.853])/365.25;
    N = length(l);
    h_tot = zeros(N, 1);
    out = exp(-( h_i(1)*sum([diff(a(l<=l_i(1)));1])+h_i(2)*sum([diff(a(l>l_i(1)&l<=l_i(2)));1])...
                +h_i(3)*sum([diff(a(l>l_i(2)&l<=l_i(3)));1])+h_i(4)*sum([diff(a(l>l_i(3)&l<=l_i(4)));1])...
                +h_i(5)*sum([diff(a(l>l_i(4)));1]) ));
end