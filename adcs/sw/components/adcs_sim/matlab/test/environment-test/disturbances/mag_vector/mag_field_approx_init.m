%% Magnetic vector unit test init file
% ----------------------------------------------------------------------- % 
% Compute data over 2 orbits to approximate magnetic field
%
% UW HuskySat-1, ADCS Subsystem
%  Last Update: T. Reynolds 3.31.18
% ----------------------------------------------------------------------- % 
run('sim_init.m')

w_prec  = 7.29211514670698e-05;
t_end   = 10800;

% Simulation parameters
run_time    = num2str(t_end);
mdl         = 'mag_field_approx';
load_system(mdl);
set_param(mdl, 'StopTime', run_time);
sim(mdl);

%% do least squares fit
A   = []; b = [];
for k = 1:length(tout)
    if( mod(k,10) == 0 )
        t = JD_ut1_J2000.Data(k);
        a = [ 1 cos(w_prec*t) sin(w_prec*t) cos(2*w_prec*t) sin(2*w_prec*t) cos(3*w_prec*t) sin(3*w_prec*t) ];
        A   = vertcat(A,kron(a,eye(3)));
        b   = vertcat(b,reshape(B_eci_T.Data(k,:),3,1));
    else
        continue;
    end
end

x   = pinv(A)*b;

% Test error with a random time
b0  = x(1:3); b1 = x(4:6); b2 = x(7:9); b3 = x(10:12); b4 = x(13:15);
b5  = x(16:18); b6  = x(19:21);

k       = round(length(tout)/2);
t       = JD_ut1_J2000.Data(k);
B_true  = B_eci_T.Data(k,:);
B_est   = @(T)(b0 + cos(w_prec.*T).*b1 + sin(w_prec.*T).*b2 + ...
            cos(2*w_prec.*T).*b3 + sin(2*w_prec.*T).*b4 + ...
            cos(3*w_prec.*T).*b5 + sin(3*w_prec.*T).*b6 );
        
err     = norm(B_true-B_est(t)');  

% Plot results over the orbits
BB  = zeros(length(B_eci_T.Time),3);
for k = 1:length(B_eci_T.Time)
    t = JD_ut1_J2000.Data(k);
    BB(k,:)     = B_est(t);
end
figure(1)
subplot(3,1,1), hold on
plot(B_eci_T.Time,B_eci_T.Data(:,1),'LineWidth',1)
plot(B_eci_T.Time,BB(:,1),'LineWidth',1)
subplot(3,1,2), hold on
plot(B_eci_T.Time,B_eci_T.Data(:,2),'LineWidth',1)
plot(B_eci_T.Time,BB(:,2),'LineWidth',1)
subplot(3,1,3), hold on
plot(B_eci_T.Time,B_eci_T.Data(:,3),'LineWidth',1)
plot(B_eci_T.Time,BB(:,3),'LineWidth',1)
xlabel('Time [s]')

if( err < 1e-7 )
    save('mag_field_approx.mat','x')
end



