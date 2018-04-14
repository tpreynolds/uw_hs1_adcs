clear variables; close all; clc;
addpath(genpath(pwd))
addpath(genpath('../../../../adcs_fsw/matlab/'))
addpath(genpath('../../../../adcs_sim/matlab/include'))
addpath(genpath('../../../../adcs_sim/matlab/tools'))
addpath(genpath('../../../../adcs_sim/matlab/lib/local/actuator'))
addpath(genpath('../../../../adcs_sim/matlab/lib/local/dynamics'))
addpath(genpath('../../../../adcs_sim/matlab/lib/local/sensor/gps'))

% Load bus stub definitions
load('bus_definitions.mat')
load('bus_definitions_fsw.mat')

% Load parameters for both flight software and simulation
fsw_params              = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);
% fsw_params.bdot         = init_bdot_controller(fsw_params);

