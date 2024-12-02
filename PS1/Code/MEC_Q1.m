% MEC
% Q1E
clear;

% Parameters
m = 1;
mu = 0.5;
k = 5;
t = 5;

% Time span
tstep = 0.1;
t_vector = 0:tstep:t;

% Initial conditions
x0 = [0; 1];

% Populate A and B matrix
A = zeros(2,2);
A(1,1) = 0;
A(1,2) = 1;
A(2,1) = -k/m;
A(2,2) = -mu/m;

B = [0 ; 1/m];

% Solve for unforced linear system
x = [];
for time = 0:tstep:t
    curr_x = expm(A*time)*x0;
    x = [x, curr_x];
end

% Plot
figure;

plot(t_vector, x(1,:));
hold on
plot(t_vector, x(2,:));
title("Plot of unforced system vs. time for spring-mass-damper system");
legend("Position of mass","Velocity of mass");
xlabel("time");
ylabel("x-value");
hold off

% Q1F
% Desired eigenvalues
p = [complex(-1,1);
    complex(-1,-1)];

% K matrix
K = place(A,B,p);
eigs = eig(A-B*K);

% Q1G
% Parameters
t = 10;
x0 = [1; 1];

% Time span
tstep = 0.1;
t_vector = 0:tstep:t;

% Solve for linear state feedback system
x = [];
for time = 0:tstep:t
    curr_x = expm((A-B*K)*time)*x0;
    x = [x, curr_x];
end

% Plot
figure;

plot(t_vector, x(1,:));
hold on
plot(t_vector, x(2,:));
title("Plot of linear state feedback system vs. time for spring-mass-damper system");
legend("Position of mass","Velocity of mass");
xlabel("time");
ylabel("x-value");
hold off