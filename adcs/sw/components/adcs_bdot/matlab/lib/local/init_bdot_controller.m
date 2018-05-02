%% BDOT Initialization File for FSW
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Subsystem
%
% Loads the parameters of the b-dot controller using predefined fsw_params.
%
% T.Reynolds 4.7.18
% ----------------------------------------------------------------------- %
function bdot = init_bdot_controller( fsw_params )

% Initial conditions
bdot.ic.RT_mt_on        = 0;
bdot.ic.RT_b_meas_valid = 0;
bdot.ic.RT_ppt_on       = 0;
bdot.ic.RT_b_body_T     = [0 0 0]';
bdot.ic.RT_dig_val      = [0 0 0]';
bdot.ic.derivative      = 0;
bdot.ic.unit_delay      = [0 0 0]';
bdot.ic.invalid_input   = [0 0 0]';

% Make sure gains are negative
bdot.gain_matrix = diag([...
            - fsw_params.actuators.magnetorquer.max_dipole_x/1.5e-6,...
            - fsw_params.actuators.magnetorquer.max_dipole_y/1.5e-6,...
            - fsw_params.actuators.magnetorquer.max_dipole_z/1.7e-6]);
                               
% NOTE: 3e-6 is the cut off value of \dot{B} below which we no longer want 
% to be saturating the torque rods.

bdot.sample_time_s  = 1/10; %[s] - sampling rate of the library

% LPF coefficients
bdot.cutoff_freq    = 2*pi*0.1; % [rad/s]
bdot.continuous_lpf = tf([bdot.cutoff_freq],[1 bdot.cutoff_freq]);
bdot.discrete_lpf   = c2d(bdot.continuous_lpf,bdot.sample_time_s);
[bdot.filter_num,bdot.filter_den]   = tfdata(bdot.discrete_lpf,'v');
% Extract second component only for use in filter
bdot.filter_num     = bdot.filter_num(2);
bdot.filter_den     = bdot.filter_den(2);

% Conversions
bdot.digital_value  = 100; 
bdot.dv_2_m_X   = fsw_params.actuators.magnetorquer.max_dipole_x/...
                    bdot.digital_value;
bdot.m_2_dv_X   = 1/bdot.dv_2_m_X;

bdot.dv_2_m_Y   = fsw_params.actuators.magnetorquer.max_dipole_y/...
                    bdot.digital_value;
bdot.m_2_dv_Y   = 1/bdot.dv_2_m_Y;

bdot.dv_2_m_Z   = fsw_params.actuators.magnetorquer.max_dipole_z/...
                    bdot.digital_value;
bdot.m_2_dv_Z   = 1/bdot.dv_2_m_Z;
 
% Tumbling status flag
bdot.B_min                  = 35e-6; % roughly min B over orbit [ T ]
bdot.omega_radps_thresh     = fsw_params.bus.omega_radps_thresh;
bdot.bdot_thresh.max        = ((sqrt(2)/2) * bdot.B_min) * ...
                                bdot.omega_radps_thresh.max;
bdot.bdot_thresh.min        = ((sqrt(2)/2) * bdot.B_min) * ...
                                bdot.omega_radps_thresh.min;
                            
end

