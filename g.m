clc; clear; close all;

% System parameters
pendulumMass = 0.1;       % Mass of pendulum bob (kg)
cartMass = 1;             % Mass of cart (kg)
pendulumLength = 0.5;     % Length of pendulum rod (m)
gravity = 9.81;           % Gravitational acceleration (m/s²)

% Linearized state-space matrices around θ = 0
A = [0, 1; 
     (cartMass + pendulumMass) * gravity / (pendulumLength * cartMass), 0];
B = [0; -1 / (pendulumLength * cartMass)];
C = [1, 0];  % Measure angular position
D = 0;       % No direct feedthrough

% Create open-loop system
openLoopSys = ss(A, B, C, D);

% PD controller design
dampingRatio = 0.5913;        % Desired damping ratio
naturalFrequency = 4.4 / dampingRatio;  % Natural frequency

% Calculate PD gains
Kp = -(naturalFrequency^2 + 21.582)/2;
Kd = -dampingRatio * naturalFrequency;

% Closed-loop system matrix (corrected formatting)
A_closed = [0, 1;
            ((cartMass + pendulumMass) * gravity + Kp) / (pendulumLength * cartMass), ...
            Kd / (pendulumLength * cartMass)];

% Create closed-loop system
closedLoopSys = ss(A_closed, B, C, D);

% Simulation parameters
timeSpan = linspace(0, 10, 500);  % Time vector for simulation
timeStep = timeSpan(2) - timeSpan(1);  % Time step

% Initial conditions (small perturbation from equilibrium)
initialAngle = pi/18;  % Approximately 10 degrees
initialAngleRate = 0;
initialState = [initialAngle; initialAngleRate];

% Simulate linearized closed-loop system
[linearResponse, time, linearStates] = initial(closedLoopSys, initialState, timeSpan);

% Simulate nonlinear system using Runge-Kutta integration
nonlinearStates = zeros(length(timeSpan), 2);
nonlinearStates(1, :) = initialState;

for step = 1:length(timeSpan)-1
    currentTime = timeSpan(step);
    currentState = nonlinearStates(step, :)';
    
    % Fourth-order Runge-Kutta method
    k1 = nonlinearDynamics(currentTime, currentState, pendulumMass, cartMass, ...
                          pendulumLength, gravity, Kp, Kd);
    k2 = nonlinearDynamics(currentTime + timeStep/2, currentState + timeStep/2 * k1, ...
                          pendulumMass, cartMass, pendulumLength, gravity, Kp, Kd);
    k3 = nonlinearDynamics(currentTime + timeStep/2, currentState + timeStep/2 * k2, ...
                          pendulumMass, cartMass, pendulumLength, gravity, Kp, Kd);
    k4 = nonlinearDynamics(currentTime + timeStep, currentState + timeStep * k3, ...
                          pendulumMass, cartMass, pendulumLength, gravity, Kp, Kd);
    
    % Update state
    nonlinearStates(step+1, :) = (currentState + (timeStep/6) * (k1 + 2*k2 + 2*k3 + k4))';
end

% Plot results
figure;

% Linear system response
subplot(2, 1, 1);
plot(time, linearStates(:, 1), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Linearized System Response with PD Control');
grid on;

% Nonlinear system response
subplot(2, 1, 2);
plot(time, nonlinearStates(:, 1), 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Nonlinear System Response with PD Control');
grid on;

% Add legend to distinguish the plots
legend('Linearized', 'Nonlinear', 'Location', 'best');

% Nonlinear pendulum dynamics with PD control
function stateDerivative = nonlinearDynamics(~, state, m, M, L, g, Kp, Kd)
    theta = state(1);
    thetaDot = state(2);
    
    % PD control law
    controlInput = -Kp * theta - Kd * thetaDot;
    
    % Intermediate calculations
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    denominator = M + m * sinTheta^2;
    
    % Nonlinear dynamics equations
    F = (sinTheta / (L * denominator)) * (-m * L * thetaDot^2 * cosTheta + (M + m) * g);
    G = (-cosTheta) / (L * denominator);
    
    thetaDoubleDot = F + G * controlInput;
    
    stateDerivative = [thetaDot; thetaDoubleDot];
end