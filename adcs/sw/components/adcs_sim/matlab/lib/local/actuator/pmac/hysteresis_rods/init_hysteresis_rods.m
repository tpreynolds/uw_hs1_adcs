%% Init Sim Hysteresis Rods
%   HuskySat-1, ADCS Subsystem
%   Last Update: B. Barzgaran - 4.13.18

%%
function hyst     = init_hysteresis_rods

% hysteresis rods properties
hyst.Hc = 1.59;
hyst.Bs = 0.73;
hyst.Br = 0.35;
hyst.V_hyst = 1;

% hysteresis switch
hyst.on     = hyst.Hc*(hyst.Bs/hyst.Br-1);
hyst.off    = hyst.Hc*(1-hyst.Bs/hyst.Br);

% rod alingment, unit vector
hyst.rod1.align = [1;0;0];
hyst.rod2.align = [0;1;0];