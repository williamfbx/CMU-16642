% MEC
% Q2
clear;

% Part c
num = [1 4 80];
denom =[2 17 158];

Ts = tf(num, denom);

zeros = zero(Ts);
poles = pole(Ts);

% Part e

step(Ts);

num_ol = [1 4 80];
denom_ol = [1 13 78];
Gs = tf(num_ol, denom_ol);
step(Gs);