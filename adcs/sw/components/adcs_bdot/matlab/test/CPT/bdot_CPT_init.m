%% BDOT Comprehensive Test
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Subsystem
%   T. Reynolds 4.11.18
%
% Initial angular velocities drawn from a uniform distribution between
% [0,11.5] deg/s or roughly [0,0.2] rad/s.
%
% Assumes sim_init.m has been run to set paths
%
% Toggle to save figures and data. 0 => no save, 1 => save.
% ----------------------------------------------------------------------- %
% Start fresh
clear variables; close all; clc
addpath(genpath('../../../matlab/'))
addpath(genpath('../../../../adcs_fsw/matlab/'))
addpath(genpath('../../../../adcs_sim/matlab/'))

set(0,'defaulttextinterpreter','latex');
figdir      = '../figs/';
datadir     = '../data/';
save_all    = 0;

%% Execute Test
k           = 1;
numTest     = 1;
W       = zeros(numTest,2); % initial/final angular velocity norm
T       = zeros(numTest,1); % time to detumble

% ----- Overrides ----- %
run_time    = 10800; % [s] -- roughly two orbits
% --------------------- %

for k = 1:numTest
    % Load bus stub definitions
    load('bus_definitions.mat')
    load('bus_definitions_fsw.mat')
    
    % Load parameters for both flight software and simulation
    fsw_params              = init_fsw_params();
    [sim_params,fsw_params] = init_sim_params(fsw_params);
    fsw_params.bdot         = init_bdot_controller(fsw_params);
    
    % Set initial tumble rate
    sim_params.dynamics.ic.rate_init    = -0.2 + (0.4)*rand(3,1);
    
    % Sim model
    mdl     = 'adcs_sim_main';
    load_system(mdl);
    set_param(mdl, 'StopTime', num2str(run_time));
    sim(mdl);
    
    % Save results of the test to an array
    W(k,1)  = norm(sim_params.dynamics.ic.rate_init);
    W(k,2)  = norm(omega_radps.Data(end,:));
    
    [val,tk]  = min(tumble.Data(2:end)); % ignore first time step
    if( val ~= 0 )
        % did not succeed in detumbling
        T(k) = -1; 
    else
        % time to detumble. +1 since we ignored first time step for tk
        T(k) = tumble.Time(tk+1); 
    end
    
end

% Save test data to csv file
