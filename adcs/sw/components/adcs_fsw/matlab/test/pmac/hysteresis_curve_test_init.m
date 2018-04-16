%% hysteresis curve test
%
% UW HuskySat-1, ADCS Subsystem

% Note: Assumes sim_init.m has been run
clc; close all;
set(0,'defaulttextinterpreter','latex');
%% Test 1

run_time    = '30';
mdl         = 'hysteresis_curve_test';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

input  = logsout.getElement('input').Values.Data
output  = logsout.getElement('output').Values.Data;

plot(input,output)
grid on