% Spinning With MTQs unit test init file
%
%
% UW HuskySat-1, ADCS Subsystem
% T. Reynolds -- 5.21.18

clear variables; close all; clc;
set(0,'defaulttextinterpreter','latex');
run('sim_init.m')
%% Test 1

% Overrides
numDay  = 0.1;
t_end   = numDay*86400;
sim_params.dynamics.ic.rate_init = [0; 0; 0];
m_x     = fsw_params.actuators.magnetorquer.max_dipole_x;
m_y     = fsw_params.actuators.magnetorquer.max_dipole_y;
m_z     = fsw_params.actuators.magnetorquer.max_dipole_z;
m_time  = 60; % seconds
% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'MTQ_spin_test';
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

