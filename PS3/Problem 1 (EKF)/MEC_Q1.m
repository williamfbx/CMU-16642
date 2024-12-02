% MEC
% Q1B
clear;

% Load data
load('calibration.mat');

% Measurement covariance
gps_diffs = [];

for index_y = 1:length(t_y)

    index_q = find(t_groundTruth == t_y(index_y));

    gps_meas = y(:,index_y);
    gt_state = q_groundTruth(1:2, index_q);
    
    gps_diffs = [gps_diffs, (gps_meas - gt_state)];

end

W = cov(gps_diffs');

% Process covariance
process_diffs = [];
T = 0.01;

for index_q = 1:(length(t_groundTruth) - 1)

    % Find state of robot at times k+1 and k
    q1_kp1 = q_groundTruth(1, index_q + 1);
    q1_k = q_groundTruth(1, index_q);
    q2_kp1 = q_groundTruth(2, index_q + 1);
    q2_k = q_groundTruth(2, index_q);
    q3_kp1 = q_groundTruth(3, index_q + 1);
    q3_k = q_groundTruth(3, index_q);

    % Find input vector at time k
    u1_k = u(1, index_q);
    u2_k = u(2, index_q);

    % Calculate process noise terms
    v1_from_q1 = ((q1_kp1 - q1_k) / (T * cos(q3_k))) - u1_k;
    v1_from_q2 = ((q2_kp1 - q2_k) / (T * sin(q3_k))) - u1_k;
    v2 = ((q3_kp1 - q3_k) / T) - u2_k;

    % Up to numerical errors, v1_from_q1 and v1_from_q2 are the same since
    % q and u are ground truth values
    v = [v1_from_q1; v2];

    process_diffs = [process_diffs, v];

end

V = cov(process_diffs');

% Q1C
clearvars -except V W

% Load data
load('kfData.mat');

% Initial parameters
T = 0.01;
q_hat = [0.355; -1.590; 0.682];
P = [25, 0, 0; 0, 25, 0; 0, 0, 0.154];

% Results storage
num_steps = length(t);
q_estimates = zeros(3, num_steps);
q_estimates(:, 1) = q_hat;

% EKF loop
for i = 1:(num_steps - 1)

    % Prediction step
    % Update mean
    q_estimates(1, i+1) = q_estimates(1, i) + T * u(1, i) * cos(q_estimates(3, i));
    q_estimates(2, i+1) = q_estimates(2, i) + T * u(1, i) * sin(q_estimates(3, i));
    q_estimates(3, i+1) = q_estimates(3, i) + T * u(2, i);

    % Update covariance
    F = [1, 0, -T * u(1, i) * sin(q_estimates(3, i)); 0, 1, T * u(1, i) * cos(q_estimates(3, i)); 0, 0, 1];
    Gamma = [T * cos(q_estimates(3, i)), 0; T * sin(q_estimates(3, i)), 0; 0, T];
    P = F * P * F' + Gamma * V * Gamma';

    % Update step (only if a new GPS measurement is received)
    if ismember(t(i+1), t_y)
        H = [1, 0, 0; 0, 1, 0];
        K = P * H' * inv(H * P * H' + W);
        y_meas = y(:, (i+1)/10);

        % The measurement equation is linear, so H = h
        q_estimates(:, i+1) = q_estimates(:, i+1) + K * (y_meas - H * q_estimates(:, i+1));
        P = (eye(3) - K * H) * P;
    end
end

% Plotting
figure;
hold on;
plot(q_groundtruth(1, :), q_groundtruth(2, :), 'b-', 'DisplayName', 'Ground Truth');
scatter(y(1, :), y(2, :), 'g', 'DisplayName', 'GPS Measurements');
plot(q_estimates(1, :), q_estimates(2, :), 'r-', 'DisplayName', 'EKF Estimate');

xlabel('x (m)');
ylabel('y (m)');
legend;
title('EKF Trajectory Estimation');
hold off;