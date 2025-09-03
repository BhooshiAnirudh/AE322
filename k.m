clc; clear; close all;

% Controller gains to analyze
proportionalGains = [-10.79, -11.5, -10];  % Different proportional gains to test
derivativeGain = 0;                        % Derivative gain (fixed at zero)

% Simulation parameters
startTime = 0;        % Start time (s)
endTime = 20;         % End time (s)
timeStep = 0.01;      % Time step (s)
timeVector = startTime:timeStep:endTime;  % Time vector

% Initial conditions
initialAngle = pi/18;  % Initial angle (rad) - approximately 10 degrees
initialState = [initialAngle; 0];  % Initial state [θ; dθ/dt]

% Create figure for response plots
figure;

% Simulate and plot system response for each Kp value
for gainIndex = 1:length(proportionalGains)
    currentKp = proportionalGains(gainIndex);
    
    % Define system dynamics with current Kp
    systemDynamics = @(t, state) [
        state(2);  % dθ/dt = angular velocity
        21.582 * state(1) - 2 * (-(currentKp * state(1)) - derivativeGain * state(2) + 1)  % d²θ/dt²
    ];
    
    % Initialize state history
    stateHistory = zeros(2, length(timeVector));
    stateHistory(:, 1) = initialState;
    
    % Runge-Kutta 4th order numerical integration
    for step = 1:length(timeVector)-1
        currentTime = timeVector(step);
        currentState = stateHistory(:, step);
        
        % Calculate RK4 coefficients
        k1 = systemDynamics(currentTime, currentState);
        k2 = systemDynamics(currentTime + timeStep/2, currentState + timeStep/2 * k1);
        k3 = systemDynamics(currentTime + timeStep/2, currentState + timeStep/2 * k2);
        k4 = systemDynamics(currentTime + timeStep, currentState + timeStep * k3);
        
        % Update state
        stateHistory(:, step+1) = currentState + timeStep/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
    
    % Extract angle from state history
    angleResponse = stateHistory(1, :);
    
    % Plot response
    subplot(length(proportionalGains), 1, gainIndex);
    plot(timeVector, angleResponse, 'r', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['System Response for K_p = ', num2str(currentKp)]);
end

% Create pole-zero maps for each Kp value
for gainIndex = 1:length(proportionalGains)
    currentKp = proportionalGains(gainIndex);
    
    % Define transfer function
    numerator = [2];  % Numerator coefficients
    denominator = [1, 2 * derivativeGain, 2 * currentKp];  % Denominator coefficients
    
    % Create transfer function
    systemTF = tf(numerator, denominator);
    
    % Create figure for pole-zero map
    figure;
    pzmap(systemTF);
    grid on;
    title(['Pole-Zero Map for K_p = ', num2str(currentKp)]);
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
end