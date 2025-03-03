function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_08_us1_5p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: United States model v1
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Bai, Y., Wagener, T., & Reed, P. (2009). A top-down framework for 
% watershed model evaluation and selection under uncertainty. Environmental
% Modelling & Software, 24(8), 901�916. 
% http://doi.org/10.1016/j.envsoft.2008.12.012


% Steps
% --- Practical ---
% 0. Handle inputs
%
% --- Model setup ---
% 1. Set out ODE
% 2. Set out constitutive functions
% 3. Determine smoothing
% 4. Determine numerical scheme

% --- Model use ---
% 5. Solve 

% --- Practical ---
% 6. Handle outputs

%% Setup
%%INPUTS
% Time step size 
delta_t = fluxInput.delta_t;

% Data
P     = fluxInput.precip./delta_t;          % [mm/delta_t] -> [mm/d]       
Ep    = fluxInput.pet./delta_t;             % [mm/delta_t] -> [mm/d]       
T     = fluxInput.temp;
t_end = length(P);

% Parameters 
% [name in documentation] = theta(order in which specified in parameter file)
alpha_ei    = theta(1);     % Fraction of intercepted rainfall [-]
m           = theta(2);     % Fraction forest [-]
smax        = theta(3);     % Maximum soil moisture [mm]
fc          = theta(4);     % Field capacity as fraction of smax [-]
alpha_ss    = theta(5);     % Subsurface routing delay [d-1]

%%INITIALISE MODEL STORES
S10         = storeInitial(1);       % Initial unsaturated storage
S20         = storeInitial(2);       % Initial saturated storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];                   % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_eusei   = zeros(1,t_end);
flux_eusveg  = zeros(1,t_end);
flux_eusbs   = zeros(1,t_end);
flux_esatveg = zeros(1,t_end);
flux_esatbs  = zeros(1,t_end);
flux_rg      = zeros(1,t_end);
flux_se      = zeros(1,t_end);
flux_qse     = zeros(1,t_end);
flux_qss     = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Unsaturated zone
% S2. Saturated zone

% EUSEI(alpha_ei,P(t)): interception as fraction of P
EUSEI = interception_3;

% EUSVEG(S1,S2,m,fc*(smax-S2),Ep(t),delta_t): transpiration from unsaturated zone
EUSVEG = evap_8;

% EUSBS(S1,S2,m,smax,Ep(t),delta_t): evaporation from unsaturated zone
EUSBS = evap_9;

% RG(P(t),S1,fc*(smax-S2)): drainage from unsaturated to saturated store
RG = saturation_1;

% SE(S1,fc*(smax-S2),delta_t): storage excess after store size change
SE = excess_1;

% ESATVEG(m,S2,smax,Ep(t),delta_t): transpiration from saturated zone
ESATVEG = evap_10;

% ESATBS(m,S2,smax,Ep(t),delta_t): evaporation from saturated zone
ESATBS = evap_5;

% QSE(RG(P,S1,fc*(smax-S2))+SE(S1,fc*(smax-S2),delta_t),S2,smax): saturation excess surface flow
QSE = saturation_1;

% QSS(alpha_ss,S2):
QSS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,1;
                                               1,1]);                       % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,1;
                                                  1,1],...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2) (P(t) - ...
                        EUSEI(alpha_ei,P(t)) - ...
                        EUSVEG(S1,S2,m,fc*(smax-S2),Ep(t),delta_t) - ...
                        EUSBS(S1,S2,m,smax,Ep(t),delta_t) - ...
                        RG(P(t),S1,fc*(smax-S2)) - ...
                        SE(S1,fc*(smax-S2),delta_t));                       % store 1 function with current flux values
    tmpf_S2 = @(S1,S2) (RG(P(t),S1,fc*(smax-S2)) + ...
                        SE(S1,fc*(smax-S2),delta_t) - ...
                        ESATVEG(m,S2,S1+S2,Ep(t),delta_t) - ...
                        ESATBS(m,S2,S1+S2,Ep(t),delta_t) - ...
                        QSE((RG(P(t),S1,fc*(smax-S2))+ SE(S1,fc*(smax-S2),delta_t)),S2,smax) - ...
                        QSS(alpha_ss,S2));                                  % store 2 function
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(eq_sys(1),eq_sys(2)),...    % system of storage equations
                        [S1old,S2old],...                                   % storage values on previous time step
                        fsolve_options);                                    % solver options
    
    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       % [tmp_sNew,tmp_resnorm,flag]
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                                    eq_sys(1),eq_sys(2)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old], ...                  % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end
    
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_eusei(t)   = EUSEI(alpha_ei,P(t));
    flux_eusveg(t)  = EUSVEG(tmp_sFlux(1),tmp_sFlux(2),m,fc*(smax-tmp_sFlux(2)),Ep(t),delta_t);
    flux_eusbs(t)   = EUSBS(tmp_sFlux(1),tmp_sFlux(2),m,smax,Ep(t),delta_t);
    flux_esatveg(t) = ESATVEG(m,tmp_sFlux(2),tmp_sFlux(1)+tmp_sFlux(2),Ep(t),delta_t);
    flux_esatbs(t)  = ESATBS(m,tmp_sFlux(2),tmp_sFlux(1)+tmp_sFlux(2),Ep(t),delta_t);
    flux_rg(t)      = RG(P(t),tmp_sFlux(1),fc*(smax-tmp_sFlux(2)));
    flux_se(t)      = SE(tmp_sFlux(1),fc*(smax-tmp_sFlux(2)),delta_t);
    flux_qse(t)     = QSE(flux_rg(t)+flux_se(t),tmp_sFlux(2),smax);
    flux_qss(t)     = QSS(alpha_ss,tmp_sFlux(2));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_eusei(t) - flux_eusveg(t) - ...
                           flux_eusbs(t) - flux_rg(t) - flux_se(t)) * delta_t;
    store_S2(t) = S2old + (flux_rg(t) + flux_se(t) - flux_esatveg(t) - ...
                           flux_esatbs(t) - flux_qse(t) - flux_qss(t)) * delta_t;    

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_eusei + flux_eusveg + flux_eusbs + ...
                         flux_esatveg + flux_esatbs) * delta_t;
    fluxOutput.Q      = (flux_qse + flux_qss) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.eusei   = flux_eusei * delta_t;
    fluxInternal.eusveg  = flux_eusveg * delta_t;
    fluxInternal.eusbs   = flux_eusbs * delta_t;
    fluxInternal.esatveg = flux_esatveg * delta_t;
    fluxInternal.esatbs  = flux_esatbs * delta_t;
    fluxInternal.rg      = flux_rg * delta_t;
    fluxInternal.se      = flux_se * delta_t;
    fluxInternal.qse     = flux_qse * delta_t;
    fluxInternal.qss     = flux_qss * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end


 




