% MEC
% Q1
clear;

% Parameters
gamma = 2;
alpha = 1;
beta = 1;
D = 1;
mu = 3;

% Populate A, B, and C matrix
A = zeros(4,4);
A(1,3) = 1;
A(2,4) = 1;
A(3,2) = 1;
A(3,3) = -3;
A(4,2) = 2;
A(4,3) = -3;

B = zeros(4,1);
B(3,1) = 1;
B(4,1) = 1;

C = [39.37 0 0 0];

% Check controllability
Q = [B A*B A*A*B A*A*A*B];
cont = rank(Q);

% Check observability
Q0 = [C; C*A; C*A*A; C*A*A*A];
obs = rank(Q0);

% Q matix for LQR
Qu = 10;

Qx = zeros(4,4);
Qx(1,1) = 30;
Qx(2,2) = 5;
Qx(3,3) = 5;
Qx(4,4) = 5;

% LQR
[Kc,S,P] = lqr(A,B,Qx,Qu);

% Eigenvalues of closed-loop system
eig_cl = eig(A-B*Kc);

% Design eigenvalues of error dynamics system
p0 = [complex(-5,1);
      complex(-5,-1);
      complex(-6,-1);
      complex(-6,1)];

K0 = place(A', C', p0)';

% Eigenvalues of error dynamics system
eig_ed = eig(A-K0*C);

% Tracking controller
Kf = -inv((C * inv(A - B * Kc) * B));

% Timespan
T = 0.01;
tspan = [0 20];
t_vector = 0:T:20;

% Initial conditions (first 4 values are actual state, and last 4 values
% are estimated state)
x0 = transpose([0, 0, 0, 0, 0.01, 0.01, -0.03, 0.01]);

% Run ode45 for original non-linear system
[t, x] = ode45(@(t, x) odefunnl(t, x, gamma, alpha, beta, D, mu, Kc, K0, Kf, A, B, C), t_vector, x0);

% Plotting
figure();
plot(t_vector,x(:,1));
hold on
plot(t_vector,x(:,5));
title("State dynamics of xc");
legend("Actual", "Estimated");
xlabel("time (sec)");
ylabel("state (m)");
hold off

figure();
plot(t_vector,x(:,2));
hold on
plot(t_vector,x(:,6));
title("State dynamics of phi");
legend("Actual", "Estimated");
xlabel("time (sec)");
ylabel("state (rad)");
hold off

figure();
plot(t_vector,x(:,3));
hold on
plot(t_vector,x(:,7));
title("State dynamics of xcdot");
legend("Actual", "Estimated");
xlabel("time (sec)");
ylabel("state (m/sec)");
hold off

figure();
plot(t_vector,x(:,4));
hold on
plot(t_vector,x(:,8));
title("State dynamics of phidot");
legend("Actual", "Estimated");
xlabel("time (sec)");
ylabel("state (rad/sec)");
hold off

% Create function for original non-linear ODE
function dxdt = odefunnl(t, x, gamma, alpha, beta, D, mu, Kc, K0, Kf, A, B, C)
    % Find feedback control from linearized ODE
    yd = 20 * square(0.02 * pi * t);

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x1hat = x(5);
    x2hat = x(6);
    x3hat = x(7);
    x4hat = x(8);
    x = [x1; x2; x3; x4];
    xhat = [x1hat; x2hat; x3hat; x4hat];

    u = -(Kc * xhat) + (Kf * yd);
    y = C*x;

    % Correction term for dxdt.hat
    corr_term = K0 * (y - C * xhat);

    dxdt = zeros(8,1);
    % Dynamics of actual system
    dxdt(1) = x3;
    dxdt(2) = x4;
    dxdt(3) = (-alpha*sin(x2)*beta*x4^2 + alpha*u - alpha*mu*x3 + cos(x2)*sin(x2)*D*beta)/(alpha*gamma - beta^2*cos(x2)^2);
    dxdt(4) = (-cos(x2)*sin(x2)*beta^2*x4^2 + u*cos(x2)*beta + sin(x2)*D*gamma - mu*x3*cos(x2)*beta)/(alpha*gamma - beta^2*cos(x2)^2);

    % Dynamics of estimated system
    dxdt(5) = x3hat + corr_term(1);
    dxdt(6) = x4hat + corr_term(2);
    dxdt(7) = (-alpha*sin(x2hat)*beta*x4hat^2 + alpha*u - alpha*mu*x3hat + cos(x2hat)*sin(x2hat)*D*beta)/(alpha*gamma - beta^2*cos(x2hat)^2) + corr_term(3);
    dxdt(8) = (-cos(x2hat)*sin(x2hat)*beta^2*x4hat^2 + u*cos(x2hat)*beta + sin(x2hat)*D*gamma - mu*x3hat*cos(x2hat)*beta)/(alpha*gamma - beta^2*cos(x2hat)^2) + corr_term(4);
end