% ----------------------------------------------------------------------- %
% UW HuskySat-1, ADCS Subsystem
%
%
%
% T. Reynolds
% ----------------------------------------------------------------------- %
%#ok<*VUNUS>

% Initial conditions
q0  = [1;0;0;0];
w0  = [0;0;0];
x   = [q0; w0];
u   = zeros(3,1);

% Magnetic field approximation parameters
b0  = fsw_params.constants.mag.b0;
b1  = fsw_params.constants.mag.b1;
b2  = fsw_params.constants.mag.b0;
b3  = fsw_params.constants.mag.b0;
b4  = fsw_params.constants.mag.b0;
b5  = fsw_params.constants.mag.b0;
b6  = fsw_params.constants.mag.b0;
w   = fsw_params.constants.mag.orbit_freq;
B   = @(t)(b0 + b1*cos(w*t) + b2*sin(w*t) + b3*cos(2*w*t) + b4*sin(2*w*t) + ...
            b5*cos(3*w*t) + b6*sin(3*w*t) ); % T
        
% Bus inertia
J   = fsw_params.bus.inertia;

% Time starts at 0s *SINCE* TLE epoch 
t0  = 0;     

% Linearize the system
[A,B1,B2] = lin_dynamics(x,u,J,B(t0));

% Set up system for solution to CARE
B       = [B1 B2];
C       = eye(6);
alpha   = 1;
R       = [ -alpha*eye(3) zeros(3); zeros(3) eye(3) ];

P   = care(A,B,C'*C,R);

K   = B1'*P; % our controller has the minus
    
disp(K)
