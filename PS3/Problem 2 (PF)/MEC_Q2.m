function M = pfTemplate()
% template and helper functions for 16-642 PS3 problem 2

rng(0); % initialize random number generator

b1 = [5,5]; % position of beacon 1
b2 = [15,5]; % position of beacon 2

% load pfData.mat
load('pfData.mat');

% initialize movie array
numSteps = length(t);
T = t(2) - t(1);
M(numSteps) = struct('cdata',[],'colormap',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         put particle filter initialization code here                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The grid is 20x10 with orientations from 0 to 2*pi
numParticles = 1000;
particles = [20 * rand(1, numParticles); 10 * rand(1, numParticles); 2 * pi * rand(1, numParticles)];

% Process noise covariance V and measurement noise covariance W
V = [1, 0; 0, 0.5];
W = [0.75, 0; 0, 0.75];

% here is some code to plot the initial scene
figure(1)
plotParticles(particles); % particle cloud plotting helper function
hold on
plot([b1(1),b2(1)],[b1(2),b2(2)],'s',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
drawRobot(q_groundTruth(:,1), 'cyan'); % robot drawing helper function
axis equal
axis([0 20 0 10])
M(1) = getframe; % capture current view as movie frame
pause
disp('hit return to continue')

% iterate particle filter in this loop
for k = 2:numSteps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              put particle filter prediction step here               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:numParticles
        % Generate process noise
        v = mvnrnd([0, 0], V);

        % Move particle
        particles(1, i) = particles(1, i) + T * (u(1, k) + v(1)) * cos(particles(3, i));
        particles(2, i) = particles(2, i) + T * (u(1, k) + v(1)) * sin(particles(3, i));
        particles(3, i) = particles(3, i) + T * (u(2, k) + v(2));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                put particle filter update step here                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % weight particles

    % Calculate expected measurement
    y_hat = zeros(2, numParticles);

    for i = 1:numParticles
        y_hat(1, i) = sqrt((particles(1, i) - b1(1))^2 + (particles(2, i) - b1(2))^2);
        y_hat(2, i) = sqrt((particles(1, i) - b2(1))^2 + (particles(2, i) - b2(2))^2);
    end

    % Calculate probability density of actual measurement
    weights = zeros(1, numParticles);

    for i = 1:numParticles
        weight1 = normpdf(y(1, k), y_hat(1, i), W(1, 1));
        weight2 = normpdf(y(2, k), y_hat(2, i), W(2, 2));
        weights(i) = weight1 * weight2;
    end

    % Normalize weights
    weights = weights / sum(weights);

    % Cumulative weight vector
    CW = cumsum(weights);

    % resample particles
    new_particles = zeros(3, numParticles);

    for i = 1:numParticles
        % Generate random number and find smallest index in CW greater than
        % number
        z = rand();
        index = find(CW > rand, 1);

        % Update particle
        new_particles(:, i) = particles(:, index);
    end

    particles = new_particles;

    % Get robot pose estimate by taking the average of the particles
    avg_particle = mean(particles, 2);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot particle cloud, robot, robot estimate, and robot trajectory here %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot beacon location, particle cloud, robot ground truth pose
    clf()
    plotParticles(particles); % particle cloud plotting helper function
    hold on
    plot([b1(1),b2(1)],[b1(2),b2(2)],'s',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[0.5,0.5,0.5]);
    drawRobot(q_groundTruth(:,k), 'cyan'); % robot drawing helper function
    axis equal
    axis([0 20 0 10])

    % Plot robot ground truth trajectory
    plot(q_groundTruth(1, 1:k), q_groundTruth(2, 1:k), 'k-', 'DisplayName', 'Ground Truth');

    % Plot robot pose estimate from particle cloud
    plot(avg_particle(1), avg_particle(2), 'r.', 'MarkerSize', 25);
    
    % capture current figure and pause
    M(k) = getframe; % capture current view as movie frame
    pause
    disp('hit return to continue')
        
end

% when you're ready, the following block of code will export the created 
% movie to an mp4 file
videoOut = VideoWriter('result.mp4','MPEG-4');
videoOut.FrameRate=5;
open(videoOut);
for k=1:numSteps
  writeVideo(videoOut,M(k));
end
close(videoOut);


% helper function to plot a particle cloud
function plotParticles(particles)
plot(particles(1, :), particles(2, :), 'go')
line_length = 0.1;
quiver(particles(1, :), particles(2, :), line_length * cos(particles(3, :)), line_length * sin(particles(3, :)))


% helper function to plot a differential drive robot
function drawRobot(pose, color)
    
% draws a SE2 robot at pose
x = pose(1);
y = pose(2);
th = pose(3);

% define robot shape
robot = [-1 .5 1 .5 -1 -1;
          1  1 0 -1  -1 1 ];
tmp = size(robot);
numPts = tmp(2);
% scale robot if desired
scale = 0.5;
robot = robot*scale;

% convert pose into SE2 matrix
H = [ cos(th)   -sin(th)  x;
      sin(th)    cos(th)  y;
      0          0        1];

% create robot in position
robotPose = H*[robot; ones(1,numPts)];

% plot robot
plot(robotPose(1,:),robotPose(2,:),'k','LineWidth',2);
rFill = fill(robotPose(1,:),robotPose(2,:), color);
alpha(rFill,.2); % make fill semi transparent

