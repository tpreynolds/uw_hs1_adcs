%% Init Sim PMAC Model
%   HuskySat-1, ADCS Subsystem
%   Last Update: B. Barzgaran - 4.12.18

%%
function pmac     = init_pmac

% Set params for each pmac library

% permanent magnet dipole
pmac.magnet.m = 1;

% hysteresis rods properties
pmac.hyst.Hc = 1.59;
pmac.hyst.Bs = 0.73;
pmac.hyst.Br = 0.35;
pmac.hyst.V_hyst = 1;

pmac.on     = pmac.hyst.Hc*(pmac.hyst.Bs/pmac.hyst.Br-1);
pmac.off    = pmac.hyst.Hc*(1-pmac.hyst.Bs/pmac.hyst.Br);

% permeability of free space
pmac.mu_0 = 1.25663706e-6;