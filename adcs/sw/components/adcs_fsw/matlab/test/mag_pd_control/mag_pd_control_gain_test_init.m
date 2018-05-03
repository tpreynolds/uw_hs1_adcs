%% MAG PD Control unit test init file
% Assumes sim_init.m has been run to set the paths
%
% Test 1: Basic test to make sure the controller reorients the bus. Used to
% choose the control gains. Uses either a random initial condition or
% identity quaternion, and commands a given slewing maneuver.

% Test 2: Recreating results from XXXXXX.

% UW HuskySat-1, ADCS Subsystem
%  Last Update: T. Reynolds 3.29.18
%% Load paths

clear variables; close all; clc;
set(0,'defaulttextinterpreter','latex');

run_test    = 2;

%% Test 1

if run_test == 1

% Start fresh
run('sim_init.m')

% Set sim time
t_end   = 40000;

% Overrides

% % Initial conditions
% sim_params.dynamics.ic.rate_init = 0.01*randn(3,1);
% temp    = randn(4,1);
% sim_params.dynamics.ic.quat_init    = temp./norm(temp);
% sim_params.dynamics.ic.quat_init    = [1; 0; 0; 0];

% Choose damping ratio and natural frequency
J  = fsw_params.bus.inertia;
z   = 1; % Critically damped
wn  = pi/2; % Small natural frequency

fsw_params.control.mag_pd_controller.p_gain = -diag([-0.03  -0.03  -0.15]);
fsw_params.control.mag_pd_controller.d_gain = -diag([-21  -21  -35]);

% fsw_params.control.pd_controller.p_gain  = wn^2.*J;
% fsw_params.control.pd_controller.d_gain  = 2*wn*z.*J;

% syms wn z
% sol = solve([wn^2.*J;2*wn*z.*J]==[-diag([-0.03  -0.03  -0.15]);-diag([-21  -21  -35])],[wn z]);

% Set commanded state
% eul_angle   = deg2rad(5);
% eul_axis    = [1; 0; 0];
% eul_axis    = eul_axis./norm(eul_axis);
% quat_cmd    = [cos(eul_angle/2); sin(eul_angle/2).*eul_axis];
% quat_cmd    = [1/2;1/2;1/2;1/2];
quat_cmd    = [1;0;0;0];
omega_cmd   = zeros(3,1);

% turn off mag noise
sim_params.sensors.magnetometer.noise = 0;

% choose dipole
fsw_params.actuators.magnetorquer.max_dipole_x  = 0.2;
fsw_params.actuators.magnetorquer.max_dipole_y  = 0.2;
fsw_params.actuators.magnetorquer.max_dipole_z  = 0.2;

% Digital value range
digital_value    = 127;

% Gains
fsw_params.control.cmd_processing.dv_2_m_X   = ...
                        fsw_params.actuators.magnetorquer.max_dipole_x /...
                            digital_value;
fsw_params.control.cmd_processing.m_2_dv_X   = 1/fsw_params.control.cmd_processing.dv_2_m_X;

fsw_params.control.cmd_processing.dv_2_m_Y   = ...
                        fsw_params.actuators.magnetorquer.max_dipole_y /...
                            digital_value;
fsw_params.control.cmd_processing.m_2_dv_Y   = 1/fsw_params.control.cmd_processing.dv_2_m_Y;

fsw_params.control.cmd_processing.dv_2_m_Z   = ...
                        fsw_params.actuators.magnetorquer.max_dipole_z /...
                            digital_value;
fsw_params.control.cmd_processing.m_2_dv_Z   = 1/fsw_params.control.cmd_processing.dv_2_m_Z;
% -----

%Control Gains
% fsw_params.control.mag_pd_controller.p_gain = 1/1000000;

% fsw_params.control.mag_pd_controller.d_gain = 1/1000000;

% -----
% Example from Torczynski, 2010
% fsw_params.bus.inertia = diag([0.037 0.036 0.006]);
% sim_params.bus.inertia = diag([0.037 0.036 0.006]);
init_quat = eul2quat(deg2rad(45*[1 1 1]))';

sim_params.dynamics.ic.quat_init = init_quat;
sim_params.dynamics.ic.rate_init = 1e-1*[0.1; 0.1; 0.1];

% fsw_params.control.mag_pd_controller.p_gain = -diag([-0.03  -0.03  -0.15]);
% fsw_params.control.mag_pd_controller.d_gain = -diag([-21  -21  -35]);

% -----

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'mag_pd_control_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %

quat        = logsout.getElement('<quaternion>').Values.Data;
omega       = logsout.getElement('<body_rates_radps>').Values.Data;
cmd_DV      = logsout.getElement('cmd_DV').Values.Data;
cmd_time    = logsout.getElement('cmd_DV').Values.Time;
real_dp     = logsout.getElement('dipole').Values.Data;
real_time   = logsout.getElement('dipole').Values.Time;
eul         = rad2deg(quat2eul(quat));
cmd_dp      = [fsw_params.control.cmd_processing.dv_2_m_X fsw_params.control.cmd_processing.dv_2_m_Y fsw_params.control.cmd_processing.dv_2_m_Z].*double(cmd_DV);

q_d         = quat_cmd; %fsw_params.bus.quat_commanded;
diff        = zeros(1,length(tout));
angle       = zeros(1,length(tout));
for i = 1:length(tout)
    q_diff  = quatmultiply(quatconj(q_d'),quat(i,:));
    diff(i) = norm( q_diff(2:4) ) ;
    angle(i) = rad2deg(2*acos(quat(i,1)));
end

% ----- End Analysis ----- %
% % Actual State Values
figure(1)
subplot(2,1,1)
plot(tout,quat)
title('Quaternion','FontSize',15)
subplot(2,1,2)
plot(tout,omega)
title('Angular Velocity [rad/s]','FontSize',15)
xlabel('Time [s]','FontSize',12)

% Commanded versus Applied Control Signals
figure(2)
subplot(2,1,1)
plot(cmd_time,cmd_dp)
title('Commanded Dipole [Nm]','FontSize',15)
subplot(2,1,2)
plot(real_time,real_dp)
title('Actual Dipole [Nm]','FontSize',15)

% Attitude Error 
figure(3), hold on
plot(tout,diff,'LineWidth',1)
plot(tout,0.02*ones(1,length(tout)),'k--')
%plot([ts ts],[0 1],'k--')
xlabel('Time [s]','FontSize',12)
title('Error')

% Angle Error
figure(5)
plot(tout,angle)
xlabel('Time [s]','FontSize',12)
ylabel('THE Euler Angle')

% Euler Angles
figure(4)
plot(tout,eul(:,1),tout,eul(:,2),tout,eul(:,3))
xlabel('Time [s]','FontSize',12)
legend('Z','Y','X')
title('Euler Angles')


%save('workspace-test-NAME.mat')

elseif run_test == 2
%% Test 2

% Start fresh
run('sim_init.m')

% Set sim time
t_end   = 40000;

% Overrides
quat_cmd    = [1;0;0;0];
omega_cmd   = zeros(3,1);

% turn off mag noise
sim_params.sensors.magnetometer.noise = 0;

% choose dipole
fsw_params.actuators.magnetorquer.max_dipole_x  = 0.2;
fsw_params.actuators.magnetorquer.max_dipole_y  = 0.2;
fsw_params.actuators.magnetorquer.max_dipole_z  = 0.2;

% Digital value range
digital_value    = 127;

% Gains
fsw_params.control.cmd_processing.dv_2_m_X   = ...
                        fsw_params.actuators.magnetorquer.max_dipole_x /...
                            digital_value;
fsw_params.control.cmd_processing.m_2_dv_X   = 1/fsw_params.control.cmd_processing.dv_2_m_X;

fsw_params.control.cmd_processing.dv_2_m_Y   = ...
                        fsw_params.actuators.magnetorquer.max_dipole_y /...
                            digital_value;
fsw_params.control.cmd_processing.m_2_dv_Y   = 1/fsw_params.control.cmd_processing.dv_2_m_Y;

fsw_params.control.cmd_processing.dv_2_m_Z   = ...
                        fsw_params.actuators.magnetorquer.max_dipole_z /...
                            digital_value;
fsw_params.control.cmd_processing.m_2_dv_Z   = 1/fsw_params.control.cmd_processing.dv_2_m_Z;
% -----

% -----
% Example from Torczynski, 2010
% fsw_params.bus.inertia = diag([0.037 0.036 0.006]);
% sim_params.bus.inertia = diag([0.037 0.036 0.006]);
init_quat = eul2quat(deg2rad(45*[1 1 1]))';

sim_params.dynamics.ic.quat_init = init_quat;
sim_params.dynamics.ic.rate_init = 1e-1*[0.1; 0.1; 0.1];

% -------------------------------------------------------


% Simulation parameters

[gain_p,gain_d] = meshgrid(0.1:0.1:2,0.1:0.1:2);
gains = [gain_p(:) gain_d(:)]; 
p = length(gains);

for i=1:p
    % Gains
    fsw_params.control.mag_pd_controller.p_gain = -gains(i,1);
    fsw_params.control.mag_pd_controller.d_gain = -gains(i,2);
    
    run_time    = num2str(t_end);
    mdl         = 'mag_pd_control_test';
    load_system(mdl);
    set_param(mdl, 'StopTime', run_time);
    sim(mdl);
    
end

% ----- Analyze Results ----- %
quat        = logsout.getElement('<quaternion>').Values.Data;
omega       = logsout.getElement('<body_rates_radps>').Values.Data;
cmd_DV      = logsout.getElement('cmd_DV').Values.Data;
cmd_time    = logsout.getElement('cmd_DV').Values.Time;
real_dp     = logsout.getElement('dipole').Values.Data;
real_time   = logsout.getElement('dipole').Values.Time;
eul         = rad2deg(quat2eul(quat));
cmd_dp      = [fsw_params.control.cmd_processing.dv_2_m_X fsw_params.control.cmd_processing.dv_2_m_Y fsw_params.control.cmd_processing.dv_2_m_Z].*double(cmd_DV);

q_d         = quat_cmd; %fsw_params.bus.quat_commanded;
diff        = zeros(1,length(tout));
angle       = zeros(1,length(tout));
for i = 1:length(tout)
    q_diff  = quatmultiply(quatconj(q_d'),quat(i,:));
    diff(i) = norm( q_diff(2:4) ) ;
    angle(i) = rad2deg(2*acos(quat(i,1)));
end

% ----- End Analysis ----- %
% % Actual State Values
figure(1)
subplot(2,1,1)
plot(tout,quat)
title('Quaternion','FontSize',15)
subplot(2,1,2)
plot(tout,omega)
title('Angular Velocity [rad/s]','FontSize',15)
xlabel('Time [s]','FontSize',12)

% Commanded versus Applied Control Signals
figure(2)
subplot(2,1,1)
plot(cmd_time,cmd_dp)
title('Commanded Dipole [Nm]','FontSize',15)
subplot(2,1,2)
plot(real_time,real_dp)
title('Actual Dipole [Nm]','FontSize',15)

% Attitude Error 
figure(3), hold on
plot(tout,diff,'LineWidth',1)
plot(tout,0.02*ones(1,length(tout)),'k--')
%plot([ts ts],[0 1],'k--')
xlabel('Time [s]','FontSize',12)
title('Error')

% Angle Error
figure(5)
plot(tout,angle)
xlabel('Time [s]','FontSize',12)
ylabel('THE Euler Angle')

% Euler Angles
figure(4)
plot(tout,eul(:,1),tout,eul(:,2),tout,eul(:,3))
xlabel('Time [s]','FontSize',12)
legend('Z','Y','X')
title('Euler Angles')


%save('workspace-test-NAME.mat')
end

