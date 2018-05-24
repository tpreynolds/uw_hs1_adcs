% ----------------------------------------------------------------------- %
%CPT_WEEK_IN_LIFE Run a long duration test of the sim
%
% UW HuskySat-1, ADCS Team
%
% Executes test 5.5 (A week in the life) for the sim. Assumes that
% sim_init.m has been run once already to set paths.
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
sim_params.CAN.COM2_on      = 1;
sim_params.CAN.pointing     = 0;
sim_params.CAN.ic.PPT_on    = 1;

% Load sim and set params
run_time    = 86400;
mdl         = 'adcs_sim_main_WIL';
load_system(mdl);
set_param(mdl,'StopTime', num2str(run_time));

