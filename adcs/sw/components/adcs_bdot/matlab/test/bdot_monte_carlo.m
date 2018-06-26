%% BDOT monte carlo Script
% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Subsystem
%   S. Rice 5/21/18
%
% Test 1:
%   Base line control gain and cut-off frequency values currently slated
%   for flight software. Run 
%
% Assumes sim_init.m has been run to set paths
%
% Toggle to save figures and data. 0 => no save, 1 => save.
save_all = 0; close all;

set(0,'defaulttextinterpreter','latex');
figdir  = './test/figs/';
datadir = './test/data/';
mu  = 398600.4418; %Standard gravitational parameter for the earth km3s-2

% %SWISSCUBE VALUES
% day_dec = 334;
% INC = 098.5033;
% RAAN = 067.1301;
% ECC = 0.0008911;
% AOP = 245.3514;
% MNA = 002.1593;
% MNM = 14.56154823;

trials = 100;
tic
incidence = zeros(trials,3);
power_after_tumble = [];
total_avg_W = [];
tumble_avg_W = [];
orbit_avg_W = [];
AOPs = [];
MNAs = [];
h = waitbar(0,'Please wait...');
for i = 1:trials
    fclose all;
    %OUR VALUES WITH RANDOM NOISE
    day_dec = 334;
    %Given data including 3 sigma
    INC =  normrnd(97.5,0.1/3);
    a = normrnd(575,20/3)+6378.1;

    %RAAN found with 10:30 LTDN which is a RAAN of 202.5 at vernal equinox,
    %then we can add 0.986 degrees/day up to DEC 1 2018
    RAAN = normrnd(095.9795,1e-1);
    ECC = normrnd(0.001,1e-5);
    AOP = 360*rand(1);
    MNA = 360*rand(1);
    MNM = sqrt(mu/a^3)*3600*24/(2*pi);  %Get the value in revs/day

    testTLE = hs1_TLE_gen(day_dec, INC, RAAN, ECC, AOP, MNA, MNM);
    fsw_params.env_estimation.orb_estimation.sgp4.orbit_tle = testTLE;
    fsw_params.bus.orbit_tle = testTLE;
    
    MNAs = [MNAs; MNA];
    AOPs = [AOPs; AOP];

    fsw_params                  = init_fsw_params();
    [sim_params, fsw_params]    = init_sim_params(fsw_params);
    fsw_params.bdot             = init_bdot_controller(fsw_params);

    % ----- Overrides ----- %
    quat = rand(4,1);
    sim_params.dynamics.ic.quat_init = quat/norm(quat);
    sim_params.dynamics.ic.rate_init = 0.3*rand(3,1);
    % --------------------- %

    % Simulation parameters
    run_time    = '10000';
    mdl         = 'bdot_simple_sim';
    load_system(mdl);
    set_param(mdl, 'StopTime', run_time);
    sim(mdl);


    tumble          = logsout.getElement('tumble').Values.Data;
    tumble_time     = logsout.getElement('tumble').Values.Time;
    
    SP_power_W = logsout.getElement('SP_power_W').Values.Data;
    SP_power_W_time = logsout.getElement('SP_power_W').Values.Time;
    
    detumble_time = find(~tumble(50:end),1);
    
    total_avg_W = [total_avg_W;mean(SP_power_W)];
    tumble_avg_W = [tumble_avg_W;mean(SP_power_W(1:detumble_time,:))];
    orbit_avg_W = [orbit_avg_W;mean(SP_power_W(detumble_time:end,:))];
    power_after_tumble = [power_after_tumble;SP_power_W(detumble_time,:)];
    
%     sc2sun_body = logsout.getElement('sc2sun_body').Values.Data;
%     sc2sun_body_time = logsout.getElement('sc2sun_body').Values.Time;
% 
%     sc_in_sun = logsout.getElement('sc_in_sun').Values.Data;
%     sc_in_sun_time = logsout.getElement('sc_in_sun').Values.Time;
%     
%     main_face_inc = sc_in_sun(end)*sc2sun_body(end,:)*[0 1 0]';
%     side_1_inc = sc_in_sun(end)*sc2sun_body(end,:)*[1 0 0]';
%     side_2_inc = sc_in_sun(end)*sc2sun_body(end,:)*[-1 0 0]';
% 
%     incidence(i,:) = [main_face_inc*(main_face_inc>0),side_1_inc*(side_1_inc>0),side_2_inc*(side_2_inc>0)];
    waitbar(i/trials)
end
toc

figure(1)
% title('Average power during tumble');
subplot(2,1,1)
histogram(sum(tumble_avg_W,2),20,'BinLimits',[1 6]);
xlabel('Average power (W)');
ylabel('trials out of 100');
title('Average power during tumble');

subplot(2,1,2)
histogram(sum(orbit_avg_W,2),20,'BinLimits',[1 6]);
xlabel('Average power (W)');
ylabel('trials out of 100');
title('Average power during full orbit post bdot');

figure(2)
histogram(sum(power_after_tumble,2),20,'BinLimits',[0 8]);
xlabel('Instantaneous Power Generation(W)');
ylabel('trials out of 100');
title('Power after tumble');
