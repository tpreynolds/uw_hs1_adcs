% Spinning With Disturbances unit test init file


% UW HuskySat-1, ADCS Subsystem
%  Last Update: T. Reynolds 8.31.17
% Load paths

clear variables; close all; clc;
set(0,'defaulttextinterpreter','latex');
run('sim_init.m')
%% Test 1

% Overrides
numDay  = 7;
t_end   = numDay*86400;
sim_params.dynamics.ic.rate_init = [0; 0; 0];
sim_params.CAN.override = 4; % force into LPM
fsw_params.sample_time_s = 1;
% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'disturbances_spin_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %
% Puts state values into a timeseries object
states      = logsout.getElement('states').Values;
sc_quat     = states.quaternion;
omega_radps = states.body_rates_radps;

% Compute RSS angular velocity
wRSS    = zeros(size(tout));
for k = 1:length(tout)
    wRSS(k) = norm(omega_radps.Data(k,:)) * fsw_params.constants.convert.rad2deg;
end

for k = 1:numDay
    wRSS_day(k) = wRSS(k*86400);
end

% ----- End Analysis ----- %
figure(1), hold on, grid on
plot(tout,omega_radps.Data,'linewidth',1)
plot(tout,fsw_params.bus.omega_radps_thresh.max*ones(size(tout)),'r--','LineWidth',1)
xlabel('Time [s]')
title('Angular velocity per axis')

figure(2), hold on, grid on
plot(tout,wRSS,'LineWidth',1)
xlabel('Time [s]')
title('RSS Tumble Rate [deg/s]')

figure(3), hold on, grid on
plot(1:numDay,wRSS_day,'ko','MarkerSize',3,'MarkerFaceColor','k')
xlabel('Time [days]')
title('RSS Tumble Rate [deg/s]')
