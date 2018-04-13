%% Init Sim PMAC Model
%   HuskySat-1, ADCS Subsystem
%   Last Update: B. Barzgaran - 4.13.18

%%
function pmac     = init_pmac

% Set params for each pmac library

% permanent magnet dipole
pmac.magnet.m = 1;
pmac.magnet.align = [0;0;1];
pmac.magnet.dipole = pmac.magnet.m*pmac.magnet.align;

% permeability of free space
pmac.mu_0 = 1.25663706e-6;

pmac.hyst = init_hysteresis_rods();