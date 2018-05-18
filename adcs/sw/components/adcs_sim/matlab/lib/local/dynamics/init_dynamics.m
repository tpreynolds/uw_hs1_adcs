%% Dynamics Init File for SIM
% ----------------------------------------------------------------------- %
% UW Husky-Sat 1, ADCS Team
%
% Initializes the dynamics block used in SIM. Quaternions are scalar first
% as per Matlab convention.
%
% T. Reynolds, 11.12.17
% ----------------------------------------------------------------------- %
function dynamics = init_dynamics()

% Quaternion FIRST
temp                    = randn(4,1);
dynamics.ic.quat_init   = temp./norm(temp);
dynamics.ic.quat_init   = [0.5 0.5 0.5 0.5]'; % override random starting condition
dynamics.ic.rate_init   = [0 0 0]';

% fsw_params.control.pd_controller.ic.torque = ...
%                             fsw_params.bus.inertia*dynamics.ic.rate_init;

end