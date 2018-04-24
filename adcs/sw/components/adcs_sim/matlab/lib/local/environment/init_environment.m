function fsw_params = init_environment(fsw_params)
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Team
%
% Initializes the whole environment block used in FSW only.
%
%   Last Edited: T. Reynolds, 8.10.17
% ----------------------------------------------------------------------- %

% Initialize all sub blocks
environment.sgp4        = init_sgp4(fsw_params);

% Update constant struct
fsw_params.constants.mag.orbit_freq = ...
                                environment.sgp4.orbit_tle(9) * ...
                                fsw_params.constants.convert.rev2rad * ...
                                fsw_params.constants.time.sec2day;

% Add sub-struct to the main fsw_params 
fsw_params.environment  = environment;

end