% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Team
%
% Set parameters to compute the position vector directly above the UW
% ground station
%
% T. Reynolds
% ----------------------------------------------------------------------- %

% Time data - Dec 1, 2018 00:00:00
gps_sec     = 518400;
gps_week    = 2029;

% UW data -- alt shouldn't matter since we only care about the unit vector
% direction here
UW_lat_lon  = fsw_params.gs_prediction.latlon;
sc_alt      = 0;

run_time    = 0;
mdl         = 'get_SEA_pos';
load_system(mdl);
set_param(mdl,'StopTime', num2str(run_time));
sim(mdl);

% Compute unit vector
r_SEA   = rSEA_ECI_m./norm(rSEA_ECI_m);

% Save vector to tool space
save('../../adcs_fsw/matlab/tools/r_SEA.mat','r_SEA');