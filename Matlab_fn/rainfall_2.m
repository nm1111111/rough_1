function [func] = rainfall_2(varargin)
%rainfall_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Rainfall based on a temperature threshold interval
% Constraints:  -
% @(Inputs):    In   - incoming precipitation flux [mm/d]
%               T    - current temperature [oC]
%               p1   - temperature threshold above which rainfall occurs [oC]
%               p2   - length of the mixed snow/rain interval [oC]
%
% WK, 08/10/2018

func = @(In,T,p1,p2) min(In,max(0,In.*(T-(p1-0.5*p2))/p2));

end

