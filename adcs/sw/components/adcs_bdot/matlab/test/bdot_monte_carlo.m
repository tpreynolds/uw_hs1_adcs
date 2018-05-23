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
power_log = [];
for i = 1:trials
    fclose all;
    %OUR VALUES WITH RANDOM NOISE
    day_dec = 334;
    %Given data including 3 sigma
    INC =  normrnd(97.5,0.1/3);
    a = normrnd(575,20/3)+6378.1;
    % a = 710

    %RAAN found with 10:30 LTDN which is a RAAN of 202.5 at vernal equinox,
    %then we can add 0.986 degrees/day up to DEC 1 2018
    RAAN = normrnd(095.9795,1e-1);
    ECC = normrnd(0.001,1e-5);
    AOP = 360*rand(1);
    MNA = 5*rand(1);
    MNM = sqrt(mu/a^3)*3600*24/(2*pi);  %Get the value in revs/day

    testTLE = hs1_TLE_gen(day_dec, INC, RAAN, ECC, AOP, MNA, MNM);
    fsw_params.env_estimation.orb_estimation.sgp4.orbit_tle = testTLE;
    fsw_params.bus.orbit_tle = testTLE;
    

    % run_test = 1;
    % ----------------------------------------------------------------------- %
    %% Test 1
    % if( run_test == 1 )
    % Initialization
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

    % Plot Results
    bdot_Tps = logsout.getElement('bdot_Tps').Values.Data;
    bdot_Tps_time = logsout.getElement('bdot_Tps').Values.Time;

    cmd_dipole_Am2 = logsout.getElement('cmd_dipole_Am2').Values.Data;
    cmd_dipole_Am2_time = logsout.getElement('cmd_dipole_Am2').Values.Time;

    body_rates_radps = logsout.getElement('body_rates_radps').Values.Data;
    body_rates_radps_time = logsout.getElement('body_rates_radps').Values.Time;

    b_meas_eci_T = logsout.getElement('b_meas_eci_T').Values.Data;
    b_meas_eci_time = logsout.getElement('b_meas_eci_T').Values.Time;

    sc2sun_body = logsout.getElement('sc2sun_body').Values.Data;
    sc2sun_body_time = logsout.getElement('sc2sun_body').Values.Time;

    sc_in_sun = logsout.getElement('sc_in_sun').Values.Data;
    sc_in_sun_time = logsout.getElement('sc_in_sun').Values.Time;

    tumble          = logsout.getElement('tumble').Values.Data;
    tumble_time     = logsout.getElement('tumble').Values.Time;
    
    SP_power_W = logsout.getElement('SP_power_W').Values.Data;
    SP_power_W_time = logsout.getElement('SP_power_W').Values.Time;
    
    power_log = [power_log;SP_power_W(end,:)];
    
%     main_face_inc = sc_in_sun(end)*sc2sun_body(end,:)*[0 1 0]';
%     side_1_inc = sc_in_sun(end)*sc2sun_body(end,:)*[1 0 0]';
%     side_2_inc = sc_in_sun(end)*sc2sun_body(end,:)*[-1 0 0]';
% 
%     incidence(i,:) = [main_face_inc*(main_face_inc>0),side_1_inc*(side_1_inc>0),side_2_inc*(side_2_inc>0)];

end
toc
% %Plotting Results
% figure(1)
% subplot(3,1,1)
% plot(cmd_dipole_Am2_time,cmd_dipole_Am2(:,1),'r')
% subplot(3,1,2)
% plot(cmd_dipole_Am2_time,cmd_dipole_Am2(:,2),'b')
% ylabel('Commanded Dipole Moment [A m2]','FontSize',12)
% subplot(3,1,3)
% plot(cmd_dipole_Am2_time,cmd_dipole_Am2(:,3),'k')
% xlabel('Time [s]','FontSize',12)
% if save_all == 1
%     SaveFigurePretty(gcf,strcat(figdir,'cmd_dipole_Am2_test1_png'));
%     saveas(gcf, strcat(figdir, 'cmd_dipole_Am2_test1'),'fig');
% end
% 
% figure(2)
% subplot(3,1,1)
% plot(body_rates_radps_time,body_rates_radps(:,1),'r','LineWidth',1)
% subplot(3,1,2)
% plot(body_rates_radps_time,body_rates_radps(:,2),'b','LineWidth',1)
% ylabel('Body Rates [rad/s]','FontSize',12)
% subplot(3,1,3)
% plot(body_rates_radps_time,body_rates_radps(:,3),'k','LineWidth',1)
% xlabel('Time [s]','FontSize',12)
% if save_all == 1
%     SaveFigurePretty(gcf,strcat(figdir,'body_rates_radps_test1_png'));
%     saveas(gcf, strcat(figdir, 'body_rates_radps_test1'),'fig');
% end
% 
% figure(3)
% subplot(2,1,1), hold on
% plot(bdot_Tps_time,bdot_Tps(:,1),'r')
% plot(bdot_Tps_time,bdot_Tps(:,2),'b')
% plot(bdot_Tps_time,bdot_Tps(:,3),'k')
% plot([0 bdot_Tps_time(end)],[fsw_params.bdot.bdot_thresh.max fsw_params.bdot.bdot_thresh.max],'r--','LineWidth',1)
% plot([0 bdot_Tps_time(end)],[fsw_params.bdot.bdot_thresh.min fsw_params.bdot.bdot_thresh.min],'r--','LineWidth',1)
% ylabel('$\dot{B}$ [T/s]','FontSize',12)
% xlabel('Time [s]','FontSize',12)
% subplot(2,1,2), hold on
% plot(tumble_time,tumble,'LineWidth',1)
% if save_all == 1
%     SaveFigurePretty(gcf,strcat(figdir,'bdot_Tps_test1_png'));
%     saveas(gcf, strcat(figdir, 'bdot_Tps_test1'),'fig');
% end
% 
% w       = 0.0693 * ones(3,1);
% B_avg   = [1.59212e-5 -6.1454e-6 4.0276e-5];
% 
% figure(4)
% subplot(3,1,1)
% plot(sc2sun_body_time,sc2sun_body(:,1),'r')
% subplot(3,1,2)
% plot(sc2sun_body_time,sc2sun_body(:,2),'b')
% subplot(3,1,3)
% plot(sc2sun_body_time,sc2sun_body(:,3),'k')
% ylabel('sc2sun body','FontSize',12)
% xlabel('Time [s]','FontSize',12)
% 
% 
% if save_all == 1
%     save(strcat(datadir,'workspace_test1.mat'),'-mat');
% end



% %% Test 2
% elseif( run_test == 2 )
% % Testing new gain and cut-off freq in the filter for better performance.
% 
% % Initialization
% clear varialbes; close all; clc
% 
% fsw_params = init_fsw_params();
% sim_params = init_sim_params(fsw_params);
% fsw_params.bdot  = init_bdot_controller(fsw_params);
% 
% % ----- Overrides ----- %
% sim_params.environment.avg_b = [1.59212e-5 -6.1454e-6 4.0276e-5]; % T
% 
% % Change gain
% fsw_params.bdot.gain_matrix    =  1.5*fsw_params.bdot.gain_matrix;
% % ----- End Overrides ----- %
% 
% % Change cut-off frequency
% fsw_params.bdot.cutoff_freq = 2*pi*0.1; % [rad/s]
% fsw_params.bdot.continuous_lpf = tf([fsw_params.bdot.cutoff_freq],[1 fsw_params.bdot.cutoff_freq]);
% fsw_params.bdot.discrete_lpf   = c2d(fsw_params.bdot.continuous_lpf,fsw_params.bdot.sample_time_s);
% [fsw_params.bdot.filter_num,fsw_params.bdot.filter_den] = tfdata(fsw_params.bdot.discrete_lpf,'v');
% % Extract second component only for use in filter
% fsw_params.bdot.filter_num     = fsw_params.bdot.filter_num(2);
% fsw_params.bdot.filter_den     = fsw_params.bdot.filter_den(2);
% 
% % --------------------- %
% 
% % Simulation parameters
% run_time    = '10000';
% mdl         = 'bdot_simple_sim';
% load_system(mdl);
% set_param(mdl, 'StopTime', run_time);
% sim(mdl);
% 
% % % Plot Results
% bdot_Tps = logsout.getElement('bdot_Tps').Values.Data;
% bdot_Tps_time = logsout.getElement('bdot_Tps').Values.Time;
% 
% cmd_dipole_Am2 = logsout.getElement('cmd_dipole_Am2').Values.Data;
% cmd_dipole_Am2_time = logsout.getElement('cmd_dipole_Am2').Values.Time;
% 
% body_rates_radps = logsout.getElement('body_rates_radps').Values.Data;
% body_rates_radps_time = logsout.getElement('body_rates_radps').Values.Time;
% 
% figure(1)
% subplot(3,1,1)
% plot(cmd_dipole_Am2_time,cmd_dipole_Am2(:,1),'r')
% subplot(3,1,2)
% plot(cmd_dipole_Am2_time,cmd_dipole_Am2(:,2),'b')
% ylabel('Commanded Dipole Moment [A m2]','FontSize',12)
% subplot(3,1,3)
% plot(cmd_dipole_Am2_time,cmd_dipole_Am2(:,3),'k')
% xlabel('Time [s]','FontSize',12)
% if save_all == 1
%     SaveFigurePretty(gcf,strcat(figdir,'cmd_dipole_Am2_test2_png'));
%     saveas(gcf, strcat(figdir, 'cmd_dipole_Am2_test2'),'fig');
% end
% 
% figure(2)
% subplot(3,1,1)
% plot(body_rates_radps_time,body_rates_radps(:,1),'r')
% subplot(3,1,2)
% plot(body_rates_radps_time,body_rates_radps(:,2),'b')
% ylabel('Body Rates [rad/s]','FontSize',12)
% subplot(3,1,3)
% plot(body_rates_radps_time,body_rates_radps(:,3),'k')
% xlabel('Time [s]','FontSize',12)
% if save_all == 1
%     SaveFigurePretty(gcf,strcat(figdir,'body_rates_radps_test2_png'));
%     saveas(gcf, strcat(figdir, 'body_rates_radps_test2'),'fig');
% end
% 
% figure(3)
% plot(bdot_Tps_time,bdot_Tps(:,1),'r')
% hold on
% plot(bdot_Tps_time,bdot_Tps(:,2),'b')
% plot(bdot_Tps_time,bdot_Tps(:,3),'k')
% ylabel('$\dot{B}$ [T/s]','FontSize',12)
% xlabel('Time [s]','FontSize',12)
% if save_all == 1
%     SaveFigurePretty(gcf,strcat(figdir,'bdot_Tps_test2_png'));
%     saveas(gcf, strcat(figdir, 'bdot_Tps_test2'),'fig');
% end
% 
% if save_all == 1
%     save(strcat(datadir,'workspace_test2.mat'),'-mat')
% end
% 
% end
