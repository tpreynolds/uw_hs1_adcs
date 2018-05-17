% function r_SEA_eci = getSeattleNADIR()
% Run simInit.m first

gps_time = [154818;
            2021];

sample_time_s = 0.1;

orbit_tle = [18;
            6848.29166666977;
            3.2923e-05;
            98.5033;
            67.1301;
            0.0008911;
            245.3514;
            2.1593;
            14.56154823];

SEA_Lat = 47;

% Set sim time
t_end   = 6000;

run_time    = num2str(t_end);
mdl         = 'Seattle_NADIR_pointing';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

Origin2SEA  = logsout.getElement('Origin2SEA').Values.Data;

R_SEA_eci = unique(Origin2SEA,'rows');
[n,~]     = size(R_SEA_eci);
j = 1;

for i = 1:n
    if norm(R_SEA_eci(i,:)) == 0
        continue
    end
    r_SEA_eci(j,:) = -R_SEA_eci(i,:)/norm(R_SEA_eci(i,:));
    j = j+1;
end

j = 1; % Vector to save
r_SEA = r_SEA_eci(j,:);
save r_SEA r_SEA