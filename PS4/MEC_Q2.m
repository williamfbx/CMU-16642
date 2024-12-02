% MEC
% Q2
clear;

A = [1, 0, 4;
     0, 1, 0;
     0, 0, 1];

B = [0.866, 0.5, 0;
     -0.5, 0.866, 0;
     0, 0, 1];

rigid_body = [-1  0  1  0 -1;
               1  1  0 -1 -1;
               1  1  1  1  1];

% A
rigid_body_a = A * B * rigid_body;

% B
rigid_body_b = B * A * rigid_body;

% C
rigid_body_c = B * rigid_body;

% D
rigid_body_d = A * B * rigid_body;

% E
rigid_body_e = B * A * rigid_body;


% Plots
plot_rigid_body(rigid_body, rigid_body_a, 'A fixed, B current');

plot_rigid_body(rigid_body, rigid_body_b, 'A fixed, B fixed');

plot_rigid_body(rigid_body, rigid_body_c, 'B fixed');

plot_rigid_body(rigid_body, rigid_body_d, 'B fixed, A fixed');

plot_rigid_body(rigid_body, rigid_body_e, 'B fixed, A current');


% Function to plot the rigid body
function plot_rigid_body(original_vertices, vertices, title_text)
    figure;
    fill(vertices(1, :), vertices(2, :), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    hold on;
    fill(original_vertices(1, :), original_vertices(2, :), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    title('SE2 Transformation');
    xlabel('X');
    ylabel('Y');
    title(title_text);
    grid on;
    axis equal;
end