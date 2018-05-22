function [ A,B1,B2 ] = lin_dynamics( x,u,J,B )
%LIN_DYNAMICS   Linearize the Satellite Dynamics
%
% Linearize the quaternion kinematics and Euler's equations about the
% operating point (x,u) assuming magnetic actuation. The input signal 'B'
% is the magnetic field vector in the body frame and determines how the
% control command influences the dynamics through
%   \tau = m x B
% This is taken into account directly, since we will assume 
%   m = K x
% for control design.
%
% UW Husky-Sat1, ADCS Subsystem
%
% T. Reynolds


q   = x(1:4);
w   = x(5:7);

dfq_dq  = 0.5*getQR([0; w]);
dfq_dw  = 0.5*getQL(q);
dfw_dq  = zeros(3,4);
dfw_dw  = J\(-skew(w)*J + skew(J*w)); 

dfq_dm  = zeros(4,3);
dfw_dm  = -inv(J)*skew(B);

A   = [ dfq_dq(2:4,2:4) dfq_dw(2:4,2:4);
        dfw_dq(1:3,2:4) dfw_dw(1:3,1:3) ];

B1  = [ dfq_dm(1:3,1:3);
        dfw_dm ];
    
B2  = [ zeros(3,3); 
        inv(J) ];    
    
end

function QR = getQR(q)
    QR = [ q(1) -q(2:4)'; q(2:4) q(1)*eye(3)-skew(q(2:4)) ];
end

function QL = getQL(q)
    QL  = [ q(1) -q(2:4)'; q(2:4) q(1)*eye(3)+skew(q(2:4)) ];
end
