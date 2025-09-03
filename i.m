clc; clear; close all;

% System parameters
pendulumMass = 0.1;       % Mass of pendulum bob (kg)
cartMass = 1;             % Mass of cart (kg)
pendulumLength = 0.5;     % Length of pendulum rod (m)
gravity = 9.81;           % Gravitational acceleration (m/s²)

% PD controller gains
proportionalGain = -155.65;  % Proportional gain
derivativeGain = -8.8;       % Derivative gain

% Define the system dynamics as an anonymous function
systemDynamics = @(t, state) [
    state(2);  % Angular velocity (dθ/dt)
    21.582 * state(1) - 2 * (-(proportionalGain * state(1)) - derivativeGain * state(2) + 1)  % Angular acceleration
];

% Simulation parameters
startTime = 0;      % Start time (s)
endTime = 20;       % End time (s)
timeStep = 0.01;    % Time step (s)
timeVector = startTime:timeStep:endTime;  % Time vector

% Initial conditions
initialAngle = pi/18;  % Initial angle (rad) - approximately 10 degrees
initialState = [initialAngle; 0];  % Initial state vector [θ; dθ/dt]

% Preallocate state matrix
stateHistory = zeros(2, length(timeVector));
stateHistory(:, 1) = initialState;

% Runge-Kutta 4th order integration
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

% Convert angle from radians to degrees
angleDegrees = stateHistory(1, :) * (180/pi);

% Calculate steady-state error
steadyStateError = -2 / (2 * proportionalGain + 21.582);
steadyStateErrorDegrees = steadyStateError * (180/pi);

% Display steady-state error
fprintf('Steady-State Error: %.4f degrees\n', steadyStateErrorDegrees);

% Plot results
figure;
plot(timeVector, angleDegrees, 'r', 'LineWidth', 1.5);
hold on;
yline(steadyStateErrorDegrees, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Steady-State Error');
grid on;
xlabel('Time (s)');
ylabel('\theta (degrees)');
title('System Response to Unit Step Disturbance');
legend('Simulated Response', 'Location', 'best');