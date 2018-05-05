%% BDOT Comprehensive Test
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Subsystem
%
% Initial angular velocities drawn from a uniform distribution between
% [0,11.5] deg/s or roughly [0,0.2] rad/s.
%
% Note: all timeseries data are pulled from logsout using a StopFcn
% callback in simulink. Can find this in:
% File --> Model Propert ies --> Model Properties --> Callbacks --> StopFcn
%
% Assumes sim_init.m has been run to set paths
%
% T. Reynolds 4.11.18
% 
% Produces simulation values that pertain to BDot, Sensor Processing and
% Estimation
% D. Tran 1 May 2018
% ----------------------------------------------------------------------- %

set(0,'defaulttextinterpreter','latex');
figdir      = '../figs/';
datadir     = '../data/';

%% Execute Test
k           = 1;
numTest     = 1;
W           = zeros(numTest,2); % initial/final angular velocity norm
T           = zeros(numTest,1); % time to detumble
DATA_fid    = strcat(datadir,'test_bdot_',char(date),'.','csv');
CPT_fid     = strcat(datadir,'bdot_CPT_',char(date),'.','csv');
SENSOR_fid  = strcat(datadir,'test_sensorproc_',char(date),'.','csv');
ESTIM_fid   = strcat(datadir,'test_estim_',char(date),'.','csv');
failed      = 0; % counter for number of failed tests ( if any )

% ----- Overrides ----- %
% One orbit is 1 hour 45 min
run_time    = 60 * 30 * 1; % [s] -- roughly two orbits
            % sec  min  hr
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
    % 1 May 2018: J. Chrisope dictates that the same rate_init
    static_vals = [0.1363 ; 0.0983 ; 0.1257];
    two_axis_zero_xz = [ 0 ; 0.4 ; 0];
    two_axis_zero_xy = [ 0 ; 0 ; 0.4 ];
    two_axis_zero_yz = [ 0.4 ; 0 ; 0];
    all_axis_zero = [ 0 ; 0 ; 0];
    sim_params.dynamics.ic.rate_init = two_axis_zero_xy;
    
    % Sim model
    mdl     = 'adcs_sim_main';
    load_system(mdl);   
    set_param(mdl, 'StopTime', num2str(run_time));
    sim(mdl);
    
    % Extract States
    states = logsout.getElement('states').Values;
    orbit_data = logsout.getElement('orbit_data').Values;
    
    % Simulation environment data
    omega_radps = states.body_rates_radps;
    env_mag_unaligned = orbit_data.mag_vec_eci;
    env_sun_unaligned = orbit_data.sun_vec_eci;
    
    % Extract bdot values for test_bdot.csv
    % Extract Actuator measurements
    act_meas = logsout.getElement('act_meas').Values;
    % Extract Commands
    cmds = logsout.getElement('commands').Values;
    MT_dv = cmds.MT_dv;
    
    % Extract sensor processing values for test_sensorproc.csv
    sens_meas = logsout.getElement('sens_meas').Values;
    sp2fsw = logsout.getElement('sp2fsw').Values;
    % Input
    mag1_vec = sens_meas.mag1_vec_body_T;
    mag2_vec = sens_meas.mag2_vec_body_T;
    imu_vec = sens_meas.omega_body_radps_gyro;
    sun_vec = sens_meas.sun_vec_body_sunsensor;
    % Output
    mag_meas = sp2fsw.mag_vec_body_T; % Also used for BDot
    gyro_meas = sp2fsw.gyro_omega_body_radps;
    sun_meas = sp2fsw.sun_vec_body_sunsensor;
    
    % Extract B-dot tumble status
    tumble = logsout.getElement('tumble').Values;
    tumble_t = logsout.getElement('tumble').Values.Time;
    
    % Extract values for test_estim.csv
    CAN_IN = logsout.getElement('CAN_IN').Values;
    envest_2_fsw = logsout.getElement('envest_2_fsw').Values;
    
    % Estim Inputs
    tle_unaligned = CAN_IN.orbit_TLE;
    met_unaligned = CAN_IN.MET;
    epoch_unaligned = CAN_IN.MET_epoch;
    
    % Estim Outputs
    sc_in_sun_unaligned = envest_2_fsw.sc_in_sun;
    sc_above_gs_unaligned = envest_2_fsw.sc_above_gs;
    sc2_sun_unit_unaligned = envest_2_fsw.sc2sun_unit;
    mag_eci_unit_unaligned = envest_2_fsw.mag_eci_unit;
    
    % Extract States
    states = logsout.getElement('states').Values;
    sc_quat = states.quaternion;
    omega_radps = states.body_rates_radps;
    % Extract Actuator measurements
    act_meas = logsout.getElement('act_meas').Values;
    % Extract Commands
    cmds = logsout.getElement('commands').Values;
    MT_dv = cmds.MT_dv;
    % Extract sensor processing output
    sp2fsw = logsout.getElement('sp2fsw').Values;
    mag_meas = sp2fsw.mag_vec_body_T;
    % Extract B-dot tumble status
    tumble = logsout.getElement('tumble').Values;
    tumble_t = logsout.getElement('tumble').Values.Time;
    
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
    i_of_sp2fsw = 2;
    i_of_CAN_IN = 2;
    i_of_envest_2_fsw = 2;
    [rows, cols] = size(tout);
    rows_tout = rows;
    mag_readings = zeros(rows, 3);
    env_mag = zeros(rows, 3);
    gyro_meas_aligned = zeros(rows, 3);
    sun_meas_aligned = zeros(rows, 3);
    % Environmental Data
    mag_data = mag_meas.Data(:,1:3);
    env_mag_data = env_mag_unaligned.Data(:,1:3);
    env_sun_data = env_sun_unaligned.Data(:,1:3);
    % Bdot Output
    mag_cmd_scaled_data = zeros(rows, 3);
    mag_volt1_data = zeros(rows, 3);
    mag_volt2_data = zeros(rows, 3);
    mag_pwm1_data = zeros(rows, 3);
    mag_pwm2_data = zeros(rows, 3);
    
    % Sensor Proc Input
    mag1_vec_data = mag1_vec.Data(:,1:4);
    mag2_vec_data = mag2_vec.Data(:,1:4);
    imu_vec_data = imu_vec.Data(:,1:4);
    sun_vec_data = sun_vec.Data(:,1:3);
    % Sensor Proc Output
    gyro_meas_data = gyro_meas.data(:,1:3);
    sun_meas_data = sun_meas.data(:,1:3);
    
    % Estim Inputs
    tle_data = tle_unaligned.Data(:,1:9);
    tle_aligned = zeros(rows, 9);
    met_data = met_unaligned.Data(:,1:1);
    met_aligned = zeros(rows, 1);
    epoch_data = epoch_unaligned.Data(:,1:1);
    epoch_aligned = zeros(rows, 1);
    % Estim Outputs
    sc_in_sun_data = sc_in_sun_unaligned.Data(:,1:1);
    sc_in_sun_aligned = zeros(rows, 1);
    sc_above_gs_data = sc_above_gs_unaligned.Data(:,1:1);
    sc_above_gs_aligned = zeros(rows, 1);
    sc2_sun_unit_data = sc2_sun_unit_unaligned.Data(:,1:3);
    sc2_sun_unit_aligned = zeros(rows, 3);
    mag_eci_unit_data = mag_eci_unit_unaligned.Data(:,1:3);
    mag_eci_unit_aligned = zeros(rows, 3);
        
    for t = 1:rows
        % Align simulation data correspondingly by time
        
        % Align sensorproc data
        if tout(t) >= sun_meas.Time(i_of_sp2fsw)
            env_mag(t,:) = env_mag_data(i_of_sp2fsw,:);
            mag_readings(t,:) = mag_data(i_of_sp2fsw,:);
            gyro_meas_aligned(t,:) = gyro_meas_data(i_of_sp2fsw,:);
            sun_meas_aligned(t,:) = sun_meas_data(i_of_sp2fsw,:);
            i_of_sp2fsw = i_of_sp2fsw + 1;
        else
            sun_meas_aligned(t,:) = sun_meas_data(i_of_sp2fsw - 1,:);
            gyro_meas_aligned(t,:) = gyro_meas_data(i_of_sp2fsw - 1,:);
            mag_readings(t,:) = mag_data(i_of_sp2fsw - 1,:);
            env_mag(t,:) = env_mag_data(i_of_sp2fsw - 1,:);
        end
        
        % Magnetometer command conversion
        % Scales
        mag_val_x = MT_dv.Data(t,:) * (100.0/127.0);
        mag_cmd_scaled_data(t,:) = mag_val_x;
        [row, cols] = size(mag_val_x);
        duty1_all_axis = zeros(1,3);
        duty2_all_axis = zeros(1,3);
        for val = 1:cols
           duty1 = 0;
           duty2 = 0;
           if (mag_val_x(val) >= 0)
               duty1 = 100 - mag_val_x(val);
           else
               duty1 = 100;
           end
           if (mag_val_x(val) < 0)
               duty2 = 100 - (-mag_val_x(val));
           else
               duty2 = 100;
           end
           duty1_all_axis(val) = duty1;
           duty2_all_axis(val) = duty2;
        end
        mag_pwm1_data(t,:) = duty1_all_axis;
        mag_pwm2_data(t,:) = duty2_all_axis;
        mag_volt1_data(t,:) = (duty1_all_axis / 100.0) * 4.5;
        mag_volt2_data(t,:) = (duty2_all_axis / 100.0) * 4.5;
        
        % Estim values alignment
        % Input alignment
        if tout(t) >= tle_unaligned.Time(i_of_CAN_IN)
            tle_aligned(t,:) = tle_data(i_of_CAN_IN,:);
            met_aligned(t,:) = met_data(i_of_CAN_IN,:);
            epoch_aligned(t,:) = epoch_data(i_of_CAN_IN,:);
            i_of_CAN_IN = i_of_CAN_IN + 1;
        else
            tle_aligned(t,:) = tle_data(i_of_CAN_IN - 1,:);
            met_aligned(t,:) = met_data(i_of_CAN_IN - 1,:);
            epoch_aligned(t,:) = epoch_data(i_of_CAN_IN - 1,:);
        end
        % Output alignment
        if tout(t) >= sc_in_sun_unaligned.Time(i_of_envest_2_fsw)
            sc_in_sun_aligned(t,:) = sc_in_sun_data(i_of_envest_2_fsw,:);
            sc_above_gs_aligned(t,:) = sc_above_gs_data(i_of_envest_2_fsw,:);
            sc2_sun_unit_aligned(t,:) = sc2_sun_unit_data(i_of_envest_2_fsw,:);
            mag_eci_unit_aligned(t,:) = mag_eci_unit_data(i_of_envest_2_fsw,:);
            i_of_envest_2_fsw = i_of_envest_2_fsw + 1;
        else
            sc_in_sun_aligned(t,:) = sc_in_sun_data(i_of_envest_2_fsw - 1,:);
            sc_above_gs_aligned(t,:) = sc_above_gs_data(i_of_envest_2_fsw - 1,:);
            sc2_sun_unit_aligned(t,:) = sc2_sun_unit_data(i_of_envest_2_fsw - 1,:);
            mag_eci_unit_aligned(t,:) = mag_eci_unit_data(i_of_envest_2_fsw - 1,:);
        end
    end
    
    header_sim = {
        'Time [s]',
        'Angular Velocity (x) [rad/s]',
        'Angular Velocity (y) [rad/s]',
        'Angular Velocity (z) [rad/s]',
        'Magnetic Sensor (x)',
        'Magnetic Sensor (y)',
        'Magnetic Sensor (z)',
        'Sun Vector (x)',
        'Sun Vector (y)',
        'Sun Vector (z)'
    };
    header_bdot = {
        'Magnetometer Reading (x) [T]',
        'Magnetometer Reading (y) [T]',
        'Magnetometer Reading (z) [T]',
        'Magnetometer Valid Bit',
        'Magnetorquer Cmds (x)',
        'Magnetorquer Cmds (y)',
        'Magnetorquer Cmds (z)',
        'Magnetorquer Cmds Scaled (x)',
        'Magnetorquer Cmds Scaled (y)',
        'Magnetorquer Cmds Scaled (z)',
        'Magnetorquer PWM 1 (x)',
        'Magnetorquer PWM 1 (y)',
        'Magnetorquer PWM 1 (z)',
        'Magnetorquer PWM 2 (x)',
        'Magnetorquer PWM 2 (y)',
        'Magnetorquer PWM 2 (z)',
        'Magnetorquer Voltage Coil 1 (x) [V]',
        'Magnetorquer Voltage Coil 1 (y) [V]',
        'Magnetorquer Voltage Coil 1 (z) [V]',
        'Magnetorquer Voltage Coil 2 (x) [V]',
        'Magnetorquer Voltage Coil 2 (y) [V]',
        'Magnetorquer Voltage Coil 2 (z) [V]'
    };
    header_sensorproc = {
        'Magnetometer 1 Input (x)',
        'Magnetometer 1 Input (y)',
        'Magnetometer 1 Input (z)',
        'Magnetometer 1 Valid Bit',
        'Magnetometer 2 Input (x)',
        'Magnetometer 2 Input (y)',
        'Magnetometer 2 Input (z)',
        'Magnetometer 2 Valid Bit',
        'IMU Input (x)',
        'IMU Input (y)',
        'IMU Input (z)',
        'IMU Input Valid Bit',
        'Sun Sensor Input (\alpha)',
        'Sun Sensor Input (\beta)',
        'Sun Sensor Input Valid Bit',
        'Magnetometer Output (x)',
        'Magnetometer Output (y)',
        'Magnetometer Output (z)',
        'IMU Output (x)',
        'IMU Output (y)',
        'IMU Output (z)',
        'Sun Sensor Output (x)',
        'Sun Sensor Output (y)',
        'Sun Sensor Output (z)'
    };
    header_estim = {
        'tle_y',
        'tle_d',
        'tle_b',
        'raan',
        'inc',
        'ecc',
        'aop',
        'mna',
        'mnm',
        'met',
        'epoch',
        'sc_in_sun',
        'sc_above_gs',
        'sc2_sun_unit (x)',
        'sc2_sun_unit (y)',
        'sc2_sun_unit (z)',
        'mag_eci_unit (x)',
        'mag_eci_unit (y)',
        'mag_eci_unit (z)'
    };
    sim_vals = [ tout omega_radps.Data env_mag env_sun_data ];
    BDot_Vals   = [ sim_vals mag1_vec_data double(MT_dv.Data) mag_cmd_scaled_data ...
        mag_pwm1_data mag_pwm2_data mag_volt1_data mag_volt2_data ];
    estim_matrix = [ sim_vals tle_aligned met_aligned epoch_aligned sc_in_sun_aligned ...
        sc_above_gs_aligned sc2_sun_unit_aligned mag_eci_unit_aligned ];
    sensorproc_matrix = [ sim_vals mag1_vec_data mag2_vec_data imu_vec_data sun_vec_data ...
        mag_readings gyro_meas_aligned sun_meas_aligned ];
    
    % Writes the header for all the csv files.
    sim_header_arr = string(header_sim)';
    bdot_header_arr = string(header_bdot)';
    sensorproc_header_arr = string(header_sensorproc)';
    estim_header_arr = string(header_estim)';
    combined_bdot_header = [ sim_header_arr bdot_header_arr ];
    combined_sensorproc_header = [ sim_header_arr sensorproc_header_arr ];
    combined_estim_header = [ sim_header_arr estim_header_arr ];
    
    % Writes the header for the bdot file
    [~, cols] = size(combined_bdot_header);
    fid = fopen(DATA_fid,'w'); 
    for i = 1:(cols - 1)
      fprintf(fid,'%s,',combined_bdot_header(i));  
    end
    fprintf(fid,'%s\n',combined_bdot_header(cols));
    fclose(fid);
    % Writes the header for the sensorproc file
    [~, cols] = size(combined_sensorproc_header);
    fid = fopen(SENSOR_fid,'w'); 
    for i = 1:(cols - 1)
      fprintf(fid,'%s,',combined_sensorproc_header(i));  
    end
    fprintf(fid,'%s\n',combined_sensorproc_header(cols));
    fclose(fid);
    % Writes the header for estim file
    [~, cols] = size(combined_estim_header);
    fid = fopen(ESTIM_fid,'w'); 
    for i = 1:(cols - 1)
      fprintf(fid,'%s,',combined_estim_header(i));  
    end
    fprintf(fid,'%s\n',combined_estim_header(cols));
    fclose(fid);
    % Writes the data to the csv files
    dlmwrite(DATA_fid,BDot_Vals,'delimiter',',','precision',9,'-append');
    dlmwrite(SENSOR_fid,sensorproc_matrix,'delimiter',',','precision',9,'-append');
    dlmwrite(ESTIM_fid,estim_matrix,'delimiter',',','precision',9,'-append');
end

% Save CPT specific data to csv file
dat = [ T W ];
csvwrite(CPT_fid,dat);
