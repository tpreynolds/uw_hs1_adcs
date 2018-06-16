% ----------------------------------------------------------------------- %
%INIT_HARDWARE_IN_LOOP Initialize the sim for HIL testing
%
% UW HuskySat-1, ADCS Team
%
% T. Reynolds
% ----------------------------------------------------------------------- %

% Start fresh
clear variables; close all; clc

% Load bus stub definitions
load('bus_definitions.mat')
load('bus_definitions_fsw.mat')

% Load parameters for both flight software and simulation
fsw_params              = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);
fsw_params.bdot         = init_bdot_controller(fsw_params);

% Test overrides
% << Set sample time for HIL here, in seconds, so 0.1 = 10Hz.
sim_params.HIL.sample_time_s = 0.1; 

% Load sim and set params
run_time    = 86400; % << Set total test time here
mdl         = 'adcs_sim_main_HIL';
load_system(mdl);
set_param(mdl,'StopTime', num2str(run_time));

