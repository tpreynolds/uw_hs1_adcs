% UW HuskySat-1, ADCS Subsystem
%   S. Rice 5/21/18
%
% Script for generating our TLE format. Will need to be updated to get a
% full TLE
% S.Rice 5/22/18
%
% Assumes sim_init.m has been run to set paths

function TLE = hs1_TLE_gen(day_dec, INC, RAAN, ECC, AOP, MNA, MNM)
mu  = 398600.4418; %  Standard gravitational parameter for the earth

%December 1, 2018
y = 18;

epoch_year          = 2000 + y;
jd_begin_of_year    = 367*epoch_year - floor((7/4)*(epoch_year + floor(10/12))) + floor(275/9) + 1721013.5;
jd_epoch_days       = day_dec + jd_begin_of_year - 2451545;

%Taken directly from SWISSCUBE
B_star      = 0.32923; % BSTAR drag term
B_star_ex   = -4;        % BSTAR drag term (exponential)
B_star      = B_star*10^(B_star_ex); 

%SWISSCUBE values
% day_dec = 334.00000000;
% INC = 098.5033;
% RAAN = 067.1301;
% ECC = 0008911;
% AOP = 245.3514;
% MNA = 002.1593
% MNM = 14.56154823

% Uncomment for saving TLE to file, but we need more information about the
% sattelite before we can really do this
% %All values that aren't specified are fixed to SWISSCUBE values
% fileID = fopen('TLE_Dump/hs1TLE.tle','w');
% fprintf(fileID,'1 35932U 09051B   %2f%11f  .00000145  00000-0  32923-4 0 00003\n',y,day_dec);
% fprintf(fileID,'2 35932 %7f %7f %7f %7f %7f %f478955',INC,RAAN,ECC,AOP,MNA,MNM);
% fclose(fileID);

TLE = [y, jd_epoch_days, B_star, INC, RAAN, ECC, AOP, MNA, MNM]';
end