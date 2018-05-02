%% State Estimation Library Test
% Test 1:
%   Provide static inertial vectors and random initial estimator states.
%   The controller should steer the estimated state the (static) commanded
%   state. Required estimator convergencence in under 60s.

% UW HuskySat-1, ADCS Subsystem
%  Primary: M. Hudoba de Badyn -- x.x.xx
%   Secondary: T. Reynolds -- 3.6.18

% Note: Assumes sim_init.m has been run
clc; close all;
set(0,'defaulttextinterpreter','latex');
%% Test 1

% seed the random initial conditions so you get the same stuff everytime
rng(155);

fsw_params.sensor_processing.gps   = init_gps_processing(fsw_params);

% Attitude conditions
quat_cmd = randn(4,1);%-[-0.1; 0.2; 0.3; 0.4];
quat_cmd = quat_cmd/norm(quat_cmd);
sim_params.dynamics.ic.quat_init        = quat_cmd';
fsw_params.estimation.ic.quat_est_init  = -[-0.1; 0.2; 0.3; 0.4];%quat_cmd';

% Angular velocity
fsw_params.estimation.ic.rt_w_body_radps    = 0*[0.01 -0.05 -0.03]';
sim_params.dynamics.ic.rate_init            = 0*[0.1 -0.05 -0.03]';

% make the noise larger (set 1 for flight values)
var_mult_mt     = 1;
var_mult_sun    = 1;
var_mult_gyro   = 1;

% Inertial measurement vectors
mag_vec_body = [4 1 1]';
mag_vec_body = mag_vec_body/norm(mag_vec_body);

sun_vec_body = [9 -3 18]';
sun_vec_body = sun_vec_body/norm(sun_vec_body);

% body frame
A = quatToMatrix(quat_cmd);
fsw_params.estimation.ic.rt_mag_eci_est = A*mag_vec_body;%.*randn(3,1);
fsw_params.estimation.ic.rt_sun_eci_est = A*sun_vec_body;%.*randn(3,1);

% reference frame
fsw_params.estimation.ic.rt_mag_body = mag_vec_body;%.*randn(3,1);
fsw_params.estimation.ic.rt_sun_body = sun_vec_body;%.*randn(3,1);

% Run the test
run_time    = '200';
mdl         = 'state_estimation_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% Plot results
h1 = figure;
plot(omega_hat.Time, omega_hat.Data(:,1),'r--'), hold on
plot(omega_hat.Time, omega_hat.Data(:,2),'b--')
plot(omega_hat.Time, omega_hat.Data(:,3),'k--')
plot(omega_true.Time, omega_true.Data(:,1),'r')
plot(omega_true.Time, omega_true.Data(:,2),'b')
plot(omega_true.Time, omega_true.Data(:,3),'k')

h2 = figure;
plot(quat_hat.Time, quat_hat.Data(:,1),'r--'), hold on
plot(quat_hat.Time, quat_hat.Data(:,2),'b--')
plot(quat_hat.Time, quat_hat.Data(:,3),'k--')
plot(quat_hat.Time, quat_hat.Data(:,4),'g--')
plot(quat_true.Time, quat_true.Data(:,1),'r')
plot(quat_true.Time, quat_true.Data(:,2),'b')
plot(quat_true.Time, quat_true.Data(:,3),'k')
plot(quat_true.Time, quat_true.Data(:,4),'g')
plot(quat_true.Time, quat_cmd(1)*ones(size(quat_true.Time)), 'k','linewidth',2)
plot(quat_true.Time, quat_cmd(2)*ones(size(quat_true.Time)), 'k','linewidth',2)
plot(quat_true.Time, quat_cmd(3)*ones(size(quat_true.Time)), 'k','linewidth',2)
plot(quat_true.Time, quat_cmd(4)*ones(size(quat_true.Time)), 'k','linewidth',2)

h3 = figure;
plot(bias_hat)
