%% Dynamics unit test init file
% Tests the dynamics block. All disturbance torques in the environment
% block are zero'd.

% UW HuskySat-1, ADCS Subsystem
%  Last Update: E. Hansen 3.31.18
%% Run sim_init first!!!

clear variables; close all; clc;
set(0,'defaulttextinterpreter','latex');

%% Tests
test = 1;

% Load parameters for both flight software and simulation
fsw_params = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);
t_end = 500;

sim_params.dynamics.ic.quat_init = [1;0;0;0];

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'solar_panel_power_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %
P_xp    = SP_power_W.signals.values(:,1);
P_xn    = SP_power_W.signals.values(:,2);
P_yp    = SP_power_W.signals.values(:,3);

quat    = sc_quat.signals.values;

figure(1), 
subplot(2,1,1), hold on, grid on
plot(sc_quat.time,quat,'LineWidth',1)
title('Quaternion')
subplot(2,1,2), hold on, grid on
plot(eul_angs.time,eul_angs.signals.values,'LineWidth',1)
legend('x','y','z')
title('Euler Angles')

figure(2), 
subplot(2,1,1), hold on, grid on
plot(SP_power_W.time,P_xp)
plot(SP_power_W.time,P_xn)
plot(SP_power_W.time,P_yp)
legend('P_x+','P_x-','P_y+')
title('Solar panel face power [W]')
subplot(2,1,2), hold on, grid on
plot(sc2sun.time,sc2sun.signals.values,'LineWidth',1)
legend('x','y','z')
title('Spacecraft to Sun Vector (ECI)')



% ----- End Analysis ----- %
