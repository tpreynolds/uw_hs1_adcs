function mag_vec = init_magnetic_field( )
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Team
%
% Initialize all parameters for the environment modeling contained in FSW
%
% T. Reynolds 6.3.18
% ----------------------------------------------------------------------- %

% Initial conditions
mag_vec.ic.ecef_2_eci       = eye(3);
mag_vec.ic.pos_eci_m        = zeros(3,1);
mag_vec.ic.time_vec_ut1     = [2018; 12; 1; 0.0; 0.0; 0.0];
mag_vec.ic.mag_vec_eci_T    = [14491; 3674; 39514] * 1e-9; % nT - total WAG
mag_vec.ic.mag_vec_eci_unit = ...
              mag_vec.ic.mag_vec_eci_T ./ norm(mag_vec.ic.mag_vec_eci_T,2);

end