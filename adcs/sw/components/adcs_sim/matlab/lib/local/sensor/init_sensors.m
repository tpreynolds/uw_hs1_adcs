function sensors = init_sensors( )
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Team
%
% Initializes all sensor parameters for the sim
%
%   Last Edited: T. Reynolds, 8.3.17
% ----------------------------------------------------------------------- %

% Initialize all sensors
%sensors.gps             = init_gps();
sensors.gyro            = init_gyroscope();
sensors.magnetometer    = init_magnetometer();
sensors.sun_sensor      = init_sun_sensor();


end