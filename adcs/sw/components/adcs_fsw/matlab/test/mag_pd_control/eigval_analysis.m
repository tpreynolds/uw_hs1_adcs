% Eigenvalue analysis of mag PD controller
% April 28, 2018
%%
q = [1;0;0;0];
w = zeros(3,1);
J = fsw_params.bus.inertia;
syms b1 b2 b3
%b = [b1 b2 b3];
b = [1.59212e-5 -6.1454e-6 4.0276e-5]; % T
Bx = [0 b(3) -b(2);
      -b(3) 0 b(1);
      b(2) -b(1) 0];

XIq     = [-q(2) -q(3) -q(4);
            q(1) -q(2)  q(3);
            q(2)  q(1) -q(4);
           -q(3)  q(4)  q(1)];

XIw     = [  0   -w(1)  w(2)  w(1);
            w(1)   0   -w(3)  w(2);
           -w(2)  w(3)   0    w(3);
           -w(1) -w(2) -w(3)   0  ];

A = [0.5*XIw(2:end,2:end) 0.5*XIq(2:end,:);
     zeros(3,3) zeros(3,3)];

B = [zeros(3);
     (J\Bx^2)/norm(b)];
 
kp1 = 0.03;
kp2 = 0.15;
kd1 = 21;
kd2 = 35;
n = 1000;
kp = logspace(-3,0,n);
kd = logspace(-2,2,n);
R = zeros(6,n);
I = zeros(6,n);

for i=1:n
    K = [kp1.*eye(3) kd(i).*eye(3)];
    Es = eig(A+B*K);
    R(:,i) = real(Es);
    I(:,i) = imag(Es);
end

figure()
hold on
for j=1:6
    plot(R(j,:),I(j,:))
end
plot(R(:,1),I(:,1),'o','MarkerFaceColor','g')
xlabel('Real')
ylabel('Imaginary')
