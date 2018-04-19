%% Simulation for PMAC system design


% UW HuskySat-1, ADCS Subsystem
%  Last Update: B. Barzgaran 4.19.18
%% Load paths
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

% ------ Overrides ------ %

% % Magnet dipole & alignment
% Ta = 8e-8; Tg = 6e-8; Tr = 1e-8; % Environmental torques, aero, gravity, radiometric
Trms = 1e-7; % from gerhardt, calculate rms for our satellite for more accuracy
Bmin = 2e-5; % minimum field strength at 600km (Tesla)
beta_max = deg2rad(10); % desired pointing accuracy
m_min = 10*Trms/(Bmin*sin(beta_max));

align_long = [0;0;1];
align_wide = [1;0;0];

m = linspace(m_min,2*m_min,10);
% align = [linspace(align_long(1),align_wide(1),10);linspace(align_long(2),align_wide(2),10);linspace(align_long(3),align_wide(3),10)];
% % alignment vector has to be unit vector:
% for i = 1:length(align(1,:))
%     align(:,i) = align(:,i)/norm(align(:,i));
% end
sim_params.actuators.pmac.magnet.align = [0;0;1];

% Hysteresis Rods
h = 35/1000; w = 55/1000; l = 70/1000;
V_max = 1/3*(h*w*l);
V_min = 1/10*V_max;
V = linspace(V_min,V_max,10);

sim_params.actuators.pmac.hyst.rod1.align = [1;0;0];
sim_params.actuators.pmac.hyst.rod2.align = [0;1;0];

% Simulation iterations
I = length(m);
J = length(align(1,:));
K = length(V);

% Simulations
t_end   = 5400;
run_time    = num2str(t_end);
model = 'pmac_sim';
load_system(model)
for ii = 1:I
    sim_params.actuators.pmac.magnet.m = m(ii);
%     for jj = 1:J
%         sim_params.actuators.pmac.magnet.align = align(:,jj);
%         sim_params.actuators.pmac.magnet.dipole = sim_params.actuators.pmac.magnet.m*sim_params.actuators.pmac.magnet.align;
%         
%         A = [align(:,jj)';zeros(2,3)];
%         Z = null(A);
%         sim_params.actuators.pmac.hyst.rod1.align = Z(:,1);
%         sim_params.actuators.pmac.hyst.rod2.align = Z(:,2);
        for kk = 1:K
            sim_params.actuators.pmac.hyst.V_hyst = V(kk);
            % Run sim
%             simOut = sim(model, 'StopTime', run_time);
            sim(model, 'StopTime', run_time);
            % Data processing
       end
%     end
end


% Evaluation




