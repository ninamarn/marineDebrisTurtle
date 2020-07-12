function out = T_ya(varargin)
    % specifies temperature
    % ...
    % user can define any function of (y, a, pars)
    global Temp
    out = 273.15 + Temp;
end