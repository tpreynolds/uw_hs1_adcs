function [CAN, fsw_params] = init_CAN( sim_params, fsw_params )
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Team
%
% Initializes all parameters for the CAN bus library used in flight
% software to simulate different signals coming from the rest of the bus.
%
%   Last Edited: T. Reynolds, 1.24.18
% ----------------------------------------------------------------------- %

% Initial conditions
CAN.ic.override                 = 0;
CAN.ic.LowPowerMode             = 0;
CAN.ic.pointing                 = 0;
CAN.ic.COM2_on                  = 0;
CAN.ic.PPT_on                   = 0;
CAN.ic.sync_pulse_period_s      = 2;
CAN.ic.PPT_on_sync_pct          = 60;
CAN.ic.MET                      = sim_params.MET.epoch;

% Update fsw.bus initial condition
CAN.ic.all              = [ CAN.ic.override;
                            CAN.ic.LowPowerMode;
                            CAN.ic.pointing;
                            CAN.ic.COM2_on;
                            CAN.ic.PPT_on];
                        
fsw_params.bus.CAN_ic   = [ CAN.ic.override;
                            CAN.ic.LowPowerMode;
                            CAN.ic.pointing;
                            CAN.ic.COM2_on;
                            CAN.ic.PPT_on ];

% sample the CAN bus at 5 Hz
CAN.sample_time_s   = 1/5;  

% Signal values
CAN.override                = 0;
CAN.LowPowerMode            = 0;
CAN.pointing                = 0;
CAN.COM2_on                 = 0;
CAN.sync_pulse_period_s     = 2;
CAN.PPT_on_sync_pct         = 60;


end

