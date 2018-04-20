%% BDOT Comprehensive Test
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Subsystem
%
% Initial angular velocities drawn from a uniform distribution between
% [0,11.5] deg/s or roughly [0,0.2] rad/s.
%
% Note: all timeseries data are pulled from logsout using a StopFcn
% callback in simulink. Can find this in:
% File --> Model Properties --> Model Properties --> Callbacks --> StopFcn
%
% Assumes sim_init.m has been run to set paths
%
% T. Reynolds 4.11.18
% ----------------------------------------------------------------------- %
set(0,'defaulttextinterpreter','latex');
figdir      = '../figs/';
datadir     = '../data/';

%% Execute Test
k           = 1;
numTest     = 1;
W           = zeros(numTest,2); % initial/final angular velocity norm
T           = zeros(numTest,1); % time to detumble
DATA_fid    = strcat(datadir,'bdot_data_',char(date),'.','csv');
CPT_fid     = strcat(datadir,'bdot_CPT',char(date),'.','csv');
failed      = 0; % counter for number of failed tests ( if any )

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
    sim_params.dynamics.ic.rate_init    = -0.2 + 0.4*rand(3,1);
    
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
        T(k)        = -1; 
        failed      = failed + 1;
        wF{failed}   = [ tout omega_radps.Data MT_dv.Data ];
        magF{failed} = [ mag_meas.Time mag_meas.Data(:,1:3) ];
        mtF{failed}  = [ MT_dv.Time MT_dv.Data ];
    else
        % time to detumble. +1 since we ignored first time step for tk
        T(k) = tumble.Time(tk+1); 
    end
       
end

% Save full test data for one run to csv file
if( T(numTest) ~= -1 )
    M   = [tout omega_radps.Data MT_dv.Data];
    N   = [mag_meas.Time mag_meas.Data(:,1:3) ];
    dlmwrite(DATA_fid,M,'delimiter',',');
    dlmwrite(DATA_fid,N,'delimiter',',','-append');
end

% Save CPT specific data to csv file
dat = [ T W ];
csvwrite(CPT_fid,dat);
