% ===== ★ ★ ★ ★ ★ =====
% pop_Caretta_caretta.m
% ©2020 Marko Jusup
% <mjusup[at]gmail.com>
% ===== ★ ★ ★ ★ ★ =====
% This script sets up the projection matrix that is used to approximate the population
% growth rate of loggerhead turtles. Specifically, we seek the solution to the Euler-Lotka
% equation 1=∫exp(-ṙa)S(a)Ḟ(a)da, where the integral—in theory—runs from 0 to ∞. For practical
% purposes, we set an upper limit at S(a)<10^(-9). Also for practical purposes, we rewrite
% the problem in discrete form 1=∑λ^(-a)S(a)Ḟ(a), where λ=exp(ṙ) and summation goes across
% all age groups for which S(a)>=10^(-9). The discrete form of the problem can further be
% rewritten in terms of a projection (i.e., Leslie) matrix in which survival between age groups
% a and a+1 is given by the ratio S(a+1)/S(a). The fecundity for age group a is just Ḟ(a). The
% population growth rate λ then becomes the dominant eigenvalue of the projection matrix. It is
% also possible generate seasonal projection matrices, which could be useful in direct simulations
% for populations that undergo important life-style changes between seasons. This is implemented
% for illustrative purposes only.
Na=150;                                                                  % number of age groups in the population
% age groups for which S(a)<10^(-9) will be eliminated later
Ns=4;                                                                      % number of seasons, e.g., set to 4 for summer, autumn, winter, spring
if ~exist('Food', 'var')                                             % set default food abundance if it has not been set
  global Food
  Food = 0.809;
end;
if ~exist('Temp', 'var')                                            % set default temperature if it has not been set
  global Temp
  Temp = 21.8;
end;
global pars bal                                                        % make parameter values accessible to subroutines
[pars, bal] = pars_Caretta_caretta();                      % load parameter values
% also balances macrochemical reactions
% the stoichiometric coefficients are returned in variable 'bal'
% these can be used to estimate, e.g., respiration or metabolic water production
[y, b] = runDEB([], 365.25*Na);                             % runs the DEB model from age a=0 to a=Na
% initial conditions are left unspecified
% the routine calculates embryonic development
% and uses that information to setup its own initial conditions
% y0 = [u_E0, u_E0-u_Eb, l_b, u_Hb];
% u_E0 scaled egg reserve
% u_Eb scaled reserve at birth
% l_b scaled length at birth
% u_Hb scaled maturity at birth
return;
Surv=zeros(Na,1);                                                  % find age group for which survival is <10^(-9)
for a = 1 : Na
    tm=365.25*(a-1);
    ix=find(min(abs(b-tm))==abs(b-tm));
    Surv(a) = get_S(b(1:ix), y(1:ix,3));
end;
Na=find(min(abs(Surv-1e-9))==abs(Surv-1e-9)); % ignore age groups with survival <10^(-9)
Smat=cell(1,Ns);                                                    % prepare seasonal projection matrices
for s = 1 : Ns
    Smat{s}=sparse(zeros(Na));
end;
for a = 1 : Na                                                          % fill seasonal projection matrices with survival and fecundities
    for s = 1 : Ns
        tm0=365.25*(a-1)+365.25*(s-1)/Ns;             % find DEB model outputs when the animal is tm0 days old (season begins)
        ix0=find(min(abs(b-tm0))==abs(b-tm0));
        tm1=365.25*(a-1)+365.25*s/Ns;                   % find DEB model outputs when the animal is tm1 days old (season ends)
        ix1=find(min(abs(b-tm1))==abs(b-tm1));
        Surv0 = get_S(b(1:ix0), y(1:ix0,3));              % survival to day tm0
        Surv1 = get_S(b(1:ix1), y(1:ix1,3));              % survival to day tm1
        Sr=Surv1/Surv0;                                             % survival through the season
        if s<Ns                                                            % only the last season advances age group a to a+1
            Smat{s}(a,a)=Sr;
        else                                                                 % reproduction too takes place in the last season
            if a<Na
                Smat{s}(a+1,a)=Sr;
            else
                Smat{s}(a,a)=Sr;
            end;
            tm0=365.25*(a-1);                                     % find DEB model outputs when the animal is tm0 days old (year begins)
            ix0=find(min(abs(b-tm0))==abs(b-tm0));
            Smat{s}(1,a)=exp(-0.69315)*sum(get_F(y(ix0:ix1,4))); % the number of females produced from the beginning of year
        end;
    end;
end;
Amat=sparse(eye(Na));                                          % annual projection matrix is just a product of seasonal matrices
for s = 1 : Ns
    Amat*=Smat{s};
end;
uH_p = pars.thr(4);                                                 % calculate age at puberty to overlay on top of population growth rates
if uH_p<max(y(:,4))
    ap = interp1(y(:,4), b, uH_p, 'pchip');
else
    ap=NaN;
end;