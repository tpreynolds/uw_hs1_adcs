%% Magnetic vector unit test init file
%
% Test 1: Run a test to compare the LS mag vec library with the WMM model
% used.
%
%
% UW HuskySat-1, ADCS Subsystem
%  Last Update: T. Reynolds 4.11.18
%% Assumes sim_init.m has been run to set the paths

% Load parameters for both flight software and simulation
fsw_params              = init_fsw_params();
[sim_params,fsw_params] = init_sim_params(fsw_params);
fsw_params.bdot          = init_bdot_controller(fsw_params);



%% Test 1

t_end = 360;

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'mag_field_approx_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

% ----- Analyze Results ----- %
B_true  = B_eci_T.Data;
B_est   = B_eci_T_est.Data;

% Compute the norm error
err     = zeros(length(tout),1);
for k = 1:length(tout)
    err(k)  = norm(B_true(k,:) - B_est(k,:));
end

% Display average norm error
fprintf('Average norm error is: %2.7f',mean(err))

%save('workspace-disturbances-test1.mat')

