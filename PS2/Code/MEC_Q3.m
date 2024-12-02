% MEC
% Q3
clear;

% Transfer function
num = [20 17];
denom = [1 9 231 400 60];
Gs = tf(num, denom);

% Create PID controller
Kp = 50;
Ki = 10;
Kd = 1;
Tf = 0;
C = pid(Kp,Ki,Kd,Tf);

% Closed-loop transfer function
Tcl = feedback(Gs*C, 1);

% Step Response
step_Tcl = step(Tcl, 50);

% Plot
step(Tcl);

S = stepinfo(Tcl);
disp(S);

% Steady-state error
sse = step_Tcl(end)-1;
disp(sse);