% MEC
% Q2B
clear;

syms alpha beta gamma D mu u x3 xcdotdot x2 x4 phidotdot
eq1 = gamma * xcdotdot - beta * phidotdot * cos(x2) + beta * x4 * x4 * sin(x2) + mu * x3 == u;
eq2 = alpha * phidotdot - beta * xcdotdot * cos(x2) - D * sin(x2) == 0;
sol = solve(eq1, eq2, xcdotdot, phidotdot);
disp(sol.xcdotdot);
disp(sol.phidotdot);
