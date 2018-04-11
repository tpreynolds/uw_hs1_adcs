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
save_all = 0;

set(0,'defaulttextinterpreter','latex');
figdir  = '../figs/';
datadir = '../data/';

%% Execute Test
k           = 1;
numTest     = 100;
W       = zeros(numTest,1); % final angular velocity norm
T       = zeros(numTest,1); % time to detumble

fsw_params          = init_fsw_params();
sim_params          = init_sim_params(fsw_params);
fsw_params.bdot     = init_bdot_controller(fsw_params);

% ----- Overrides ----- %
sim_params.dynamics.ic.rate_init    = -0.2 + (0.4)*rand(3,1);
run_time    = 9500; % [s] -- roughly two orbits
% --------------------- %

% Simulation parameters
mdl     = 'adcs_sim_main';
load_system(mdl); 
set_param(mdl, 'StopTime', num2str(run_time)); 
sim(mdl);

% Save results of the test to an array
W(k)    = norm(omega_radps(end,:));
