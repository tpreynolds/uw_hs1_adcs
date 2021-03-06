%% Sunsensor sensor processing unit test init file

% Test 1: 

% UW HuskySat-1, ADCS Subsystem
%  Last Update: T. Reynolds 12.12.17
%% Load paths

% Start fresh
clear variables; close all; clc
set(0,'defaulttextinterpreter','latex');
% addpath(genpath('../../../../matlab/')) % adds the fsw libs
% addpath(genpath('../../../../../adcs_sim/matlab/')) % add the sim libs

% Load bus stub definitions
load('bus_definitions.mat')

run_test    = 1;
%% Test 1

if run_test == 1

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

% Overrides
t_end   = 100;
% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'sunsensor_processing_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %



% ----- End Analysis ----- %

figure(1)
subplot(2,2,1)
plot(tout(11:end),gps_time_fsw(11:end,1),'LineWidth',1)
title('FSW GPS Seconds','FontSize',15)
subplot(2,2,2)
plot(tout(11:end),gps_time_fsw(11:end,2),'LineWidth',1)
title('FSW GPS Weeks','FontSize',15)
subplot(2,2,3)
plot(gps_sensor_time(2:end),gps_time_sensor(2:end,1),'LineWidth',1)
title('Sensor GPS Seconds','FontSize',15)
xlabel('Time [s]','FontSize',12)
subplot(2,2,4)
plot(gps_sensor_time(2:end),gps_time_sensor(2:end,2),'LineWidth',1)
title('Sensor GPS Weeks','FontSize',15)
xlabel('Time [s]','FontSize',12)

%save('workspace-test1-GPS.mat')

elseif run_test == 2
%% Test 2
clear variables; close all; clc

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);

% Overrides
t_end   = 100;
% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'GPS_processing_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze + Plot Results ----- %


% ----- End Analysis ----- %

%save('workspace-test2-NAME.mat')



end

