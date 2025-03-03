function [func] = baseflow_1(~)
%baseflow_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Outflow from a linear reservoir
% Constraints:  -
% @(Inputs):    p1   - time scale parameter [d-1]
%               S    - current storage [mm]
%
% WK, 05/10/2018

func = @(p1,S) p1.*S;

end