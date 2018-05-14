%% MAG PD Control initial condition test init file
% Assumes sim_init.m has been run to set the paths

% Test 1: Test to validate initial condition stability

% UW HuskySat-1, ADCS Subsystem
%  Last Update: T. Reynolds 3.29.18
%% Load paths

clear variables; close all; clc;
set(0,'defaulttextinterpreter','latex');
%% Test 1

% Start fresh
run('sim_init.m')

% Set sim time
t_end   = 20000;

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
% init_quat = eul2quat(deg2rad(45*[1 1 1]))';
% 
% sim_params.dynamics.ic.quat_init = init_quat;
% sim_params.dynamics.ic.rate_init = 1e-1*[0.1; 0.1; 0.1];
IC_quat_eul = linspace(-60,60,10);
IC_quat     = deg2rad(IC_quat_eul);
IC_rate     = linspace(-0.2,0.2,10);
[quat_plot,rate_plot] = meshgrid(IC_quat,IC_rate);

% -------------------------------------------------------

% Simulation parameters
initcond = [quat_plot(:) rate_plot(:)]; 
p = length(initcond);
settle_t = zeros(p,3);
steady = zeros(p,1);


% Gains
p_gain = ;
d_gain = ;
fsw_params.control.mag_pd_controller.p_gain = -p_gain*eye(3);
fsw_params.control.mag_pd_controller.d_gain = -d_gain*eye(3);


for i=1:p
    % Initial conditions
    
    init_quat = eul2quat(initcond(i,1)*[1 1 1])';
    sim_params.dynamics.ic.quat_init = init_quat;
    sim_params.dynamics.ic.rate_init = initcond(i,2)*[1; 1; 1];
    
    run_time    = num2str(t_end);
    mdl         = 'mag_pd_control_test';
    load_system(mdl);
    set_param(mdl, 'StopTime', run_time);
    sim(mdl);
    avg_eul     = logsout.getElement('avg_eul_proc').Values.Data;
    avg_eul_t   = logsout.getElement('avg_eul_proc').Values.Time;
    N = length(avg_eul_t);
    settle_t(i,2) = N;
    settle_t(i,3) = avg_eul_t(end);
    for j=1:N
        if isnan(avg_eul(N-j-1)) ~= 1
            settle_t(i,1) = avg_eul_t(N-j-1);
            settle_t(i,2) = N-j-1;
        else
            break
        end
    end
    if settle_t(i,2) ~= N
        steady(i) = rms(avg_eul(settle_t(i,2):end));
    end
    i
end

ST = reshape(settle_t(:,1),size(gain_p)); %settling time
SS = reshape(steady,size(gain_p)); %steady state rms

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
figure(4)
plot(tout,angle)
xlabel('Time [s]','FontSize',12)
ylabel('THE Euler Angle')

% Euler Angles
figure(5)
plot(tout,eul(:,1),tout,eul(:,2),tout,eul(:,3))
xlabel('Time [s]','FontSize',12)
legend('Z','Y','X')
title('Euler Angles')
%%
% Settling time
figure(6)
surf(gain_p,gain_d,ST)
view([142.5,30])
xlabel('p Gain')
ylabel('d Gain')
zlabel('Settling Time')
grid minor

% Steady State rms
figure(7)
surf(gain_p,gain_d,SS)
view([142.5,30])
xlabel('p Gain')
ylabel('d Gain')
zlabel('Steady State RMS')
grid minor

%save('workspace-test-NAME.mat')
