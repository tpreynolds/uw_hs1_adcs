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
align = [linspace(align_long(1),align_wide(1),10);linspace(align_long(2),align_wide(2),10);linspace(align_long(3),align_wide(3),10)];
% alignment vector has to be unit vector:
for i = 1:length(align(1,:))
    align(:,i) = align(:,i)/norm(align(:,i));
end

% Hysteresis Rods
sim_params.actuators.pmac.hyst

% Simulation iterations
I = length(m);
J = length(align(1,:));


% Simulations
model = 'pmac_sim';
load_system(model)
for ii = 1:I
    sim_params.actuators.pmac.magnet.m = m(ii);
    for jj = 1:J
        sim_params.actuators.pmac.magnet.align = align(:,jj);
        sim_params.actuators.pmac.magnet.dipole = m*align;
        
        
        
        
        % Run sim
       simOut = sim(model); 
    end
end


% Evaluation




