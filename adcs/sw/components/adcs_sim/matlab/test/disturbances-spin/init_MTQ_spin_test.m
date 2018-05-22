% Spinning With MTQs unit test init file
%
%
% UW HuskySat-1, ADCS Subsystem
% T. Reynolds -- 5.21.18

clear variables; close all; clc;
set(0,'defaulttextinterpreter','latex');
run('sim_init.m')
%% Test 1

% Test parameters
test    = 0; % Set 0 to use pegged dipoles, 1 to use reverse b-dot
numDay  = 0.1;
t_end   = 5400;%numDay*l86400;
m_x     = fsw_params.actuators.magnetorquer.max_dipole_x;
m_y     = fsw_params.actuators.magnetorquer.max_dipole_y;
m_z     = fsw_params.actuators.magnetorquer.max_dipole_z;
m_mat   = diag([m_x,m_y,m_z]);
m_time  = 5400; % pegged actuation time

% Overrides
sim_params.dynamics.ic.rate_init    = [0; 0; 0];
fsw_params.bdot.gain_matrix         = - fsw_params.bdot.gain_matrix;
% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'MTQ_spin_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl); 
% -----

% Puts state values into a timeseries object
states      = logsout.getElement('states').Values;
sc_quat     = states.quaternion;
omega_radps = states.body_rates_radps;

% Plot Results
figure(1), hold on, grid on
plot(tout,omega_radps.Data,'linewidth',1)
plot(tout,fsw_params.bus.omega_radps_thresh.max*ones(size(tout)),'r--','LineWidth',1)
xlabel('Time [s]')
title('Angular velocity per axis')

figure(2), hold on, grid on
plot(wRSS.time,wRSS.signals.values,'LineWidth',1)
xlabel('Time [s]','FontSize',14)
title('RSS Tumble Rate [deg/s]','FontSize',16)
if( test == 1 )
axes('position',[.65 .175 .25 .25])
box on
ind = ( 0 < wRSS.time ) & ( wRSS.time < 600 );
plot(wRSS.time(ind),wRSS.signals.values(ind),'LineWidth',1)
axis tight, grid on
end


% end test
