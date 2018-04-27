%% Sunsensor sensor processing unit test init file

% Test 1: 

% UW HuskySat-1, ADCS Subsystem
% Last Update: M. Hudoba de Badyn April 20 2018
%% Load paths

% Start fresh
clear variables; close all; clc
set(0,'defaulttextinterpreter','latex');
addpath(genpath('../../../../matlab/')) % adds the fsw libs
addpath(genpath('../../../../../adcs_sim/matlab/')) % add the sim libs

% Load bus stub definitions
load('bus_definitions.mat')

fprintf('Sign tests for alpha and beta, and the x/y components of sun vector\n')

%% Test 1

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = pi/6;
beta = pi/6;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,1:3);
fprintf('Alpha and Beta are both >0\n')
fprintf('%i should be %i\n',sign(sun_vec(1)), sign(alpha))
fprintf('%i should be %i\n',sign(sun_vec(2)), sign(beta))

% ----- End Analysis ----- %


%% Test 2

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = -pi/6;
beta = pi/6;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,1:3);
fprintf('Alpha <0, Beta >0\n')
fprintf('%i should be %i\n',sign(sun_vec(1)), sign(alpha))
fprintf('%i should be %i\n',sign(sun_vec(2)), sign(beta))

% ----- End Analysis ----- %

%% Test 3

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = pi/6;
beta = -pi/6;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,1:3);
fprintf('Alpha >0, Beta <0\n')
fprintf('%i should be %i\n',sign(sun_vec(1)), sign(alpha))
fprintf('%i should be %i\n',sign(sun_vec(2)), sign(beta))

% ----- End Analysis ----- %

%% Test 4

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = -pi/6;
beta = -pi/6;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,1:3);
fprintf('Alpha and Beta are both <0\n')
fprintf('%i should be %i\n',sign(sun_vec(1)), sign(alpha))
fprintf('%i should be %i\n',sign(sun_vec(2)), sign(beta))

% ----- End Analysis ----- %


fprintf('Valid tests for alpha and beta being out of range\n')

%% Test 5

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = -pi/2;
beta = -pi/6;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,4);
fprintf('Alpha < -60deg\n')
fprintf('%i should be 0\n', sun_vec )


% ----- End Analysis ----- %

%% Test 6

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = pi/2;
beta = -pi/6;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,4);
fprintf('Alpha > 60deg\n')
fprintf('%i should be 0\n', sun_vec )


% ----- End Analysis ----- %

%% Test 7

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = -pi/6;
beta = -pi/2;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,4);
fprintf('beta < -60deg\n')
fprintf('%i should be 0\n', sun_vec )

% ----- End Analysis ----- %

%% Test 8

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

alpha = -pi/6;
beta = pi/2;
sun_vec_valid = 1;

% Overrides
t_end   = 1;

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_valid_flag_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

sun_vec    = logsout.getElement('sunsensor_vec_fsw').Values.Data;
sun_vec = sun_vec(1,4);
fprintf(' beta > 60deg\n ')
fprintf('%i should be 0\n', sun_vec )

% ----- End Analysis ----- %


