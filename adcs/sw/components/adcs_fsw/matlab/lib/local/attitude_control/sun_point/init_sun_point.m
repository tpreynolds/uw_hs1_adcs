function sun_point = init_sun_point(fsw_params)
% ----------------------------------------------------------------------- %
%INIT_SUN_POINT     Initialize the sun pointing controller
%
% UW HuskySat-1, ADCS Team
% 
% Sets all parameters for the sun pointing controller used to acquire the
% sun when the sunsensor is invalid.
% 
% T. Reynolds -- 5.9.18
% ----------------------------------------------------------------------- %

% Sample times
sun_point.sample_time_s     = fsw_params.sample_time_s;

% Initial conditions
sun_point.ic.meas_ss_valid  = 0;
sun_point.ic.meas_ss_body_unit      = [0 0 0];
sun_point.ic.meas_pd_body_unit      = [0 0 0];
sun_point.ic.meas_body_rate_radps   = [0 0 0];
sun_point.ic.error_vec      = [0 0 0];
sun_point.ic.cmd_torque_nm  = [0 0 0];
sun_point.ic.sun_body_unit  = [0;0;1];

% % Set gains
% J   = fsw_params.bus.inertia;
% z   = 1;
% wn  = 0.005 * 2 * pi;
% sun_point.prop_gain = wn^2 * J;
% sun_point.drv_gain  = 2 * wn * z * J;

sun_point.prop_gain     = ...
         [      -7.15925376030017e-16      2.95183064103284e-05      1.16471856488024e-05;
      -2.9518306409416e-05      6.91389624109808e-16       6.3973654111362e-06;
     -1.16471856520801e-05     -6.39736541290793e-06      1.52865094782168e-17];
sun_point.drv_gain = ...
         [       -1.13493695706464e-09      3.00242063478115e-05      1.16904253272767e-05;
      -3.0013255295274e-05      7.31311204065323e-10      6.42213829448181e-06;
     -1.18422341375943e-05     -6.50671819922309e-06      4.03625743741717e-10];
end