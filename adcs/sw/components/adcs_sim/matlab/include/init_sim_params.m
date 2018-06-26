% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Team
%
% Define all parameters to be used by SIM here. This is the second file to
% be called by 'SimInit.m' to initialize simulation data for the full 'Main
% Simulation.slx'. 
%   Takes in already defined fsw_params so that sim values can be defined
%   using the fsw values
%
% T.Reynolds 4.18.18
% ----------------------------------------------------------------------- %
function [sim_params,fsw_params] = init_sim_params(fsw_params)

% ----- Spacecraft Parameters ----- %
sim_params.bus  = fsw_params.bus;
sim_params.bus.decimation.long          = 60;
sim_params.bus.decimation.short         = 6;
sim_params.bus.decimation.num_pts_long  = 60000;
sim_params.bus.decimation.num_pts_short = 6000;
[sim_params.MET, fsw_params]  = init_MET(fsw_params);
% --------------------------------- %

% ----- Physical Bus Signal Emulators ----- %
[sim_params.CAN, fsw_params]  = init_CAN(sim_params, fsw_params);
% ----------------------------------------- %

% ----- Dynamics ----- %
sim_params.dynamics = init_dynamics();
% -------------------- %

% ----- Sensors ----- %
sim_params.sensors  = init_sensors();
% ------------------- %

% ----- Actuators ----- %
sim_params.actuators    = init_actuators(fsw_params);
% --------------------- %

% ----- Environment ----- %
sim_params  = init_environment_sim(fsw_params,sim_params);
% ----------------------- %

% ----- Estimation ----- %
fsw_params.estimation   = init_extended_kalman_filter(sim_params);
% ---------------------- %

 end
