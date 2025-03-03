function [func] = percolation_4(~)
%percolation_4 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Demand-based percolation scaled by available moisture
% Constraints:  f <= S/dt
%               f >= 0          prevents erratic numerical behaviour
% @(Inputs):    p1   - base percolation rate [mm/d]
%               p2   - percolation rate increase due moisture deficiencies [mm/d]
%               p3   - non-linearity parameter [-]
%               p4   - summed deficiency across all model stores [mm]
%               p5   - summed capacity of model stores [mm]
%               S    - current storage in the supplying store [mm]
%               Smax - maximum storage in the supplying store [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

% Note: for certain extreme parameter values (very small stores, highly
% non-linear p3) and small computational errors that lead to small negative
% S values, this function behaves erratically. The max(0,S/Smax) part
% prevents this. Similarly, the first max(0,...) part prevents negative
% percolation demands as a result of small numerical errors.

func = @(p1,p2,p3,p4,p5,S,Smax,dt) max(0,min(S/dt,max(0,S./Smax).*(p1.*(1+p2.*(p4./p5).^(1+p3)))));

end

