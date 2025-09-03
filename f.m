clc; clear; close all;

% System parameters
pendulumMass = 0.1;       % Mass of pendulum bob (kg)
cartMass = 1;             % Mass of cart (kg)
pendulumLength = 0.5;     % Length of pendulum rod (m)
gravity = 9.81;           % Gravitational acceleration (m/s²)

% Linearized state matrices for two equilibrium positions
A_matrices = {
    [0, 1; (cartMass + pendulumMass) * gravity / (pendulumLength * cartMass), 0],  % Around θ = 0
    [0, 1; -(cartMass + pendulumMass) * gravity / (pendulumLength * cartMass), 0]  % Around θ = π
};

B_matrices = {
    [0; -1 / (pendulumLength * cartMass)],  % Input matrix for θ = 0
    [0; 1 / (pendulumLength * cartMass)]    % Input matrix for θ = π
};

% Output matrix (measure angle only)
C = [1, 0];
D = 0;  % No direct feedthrough

% Simulation parameters
timeSpan = linspace(0, 20, 500);  % Time vector for simulation
timeStep = timeSpan(2) - timeSpan(1);  % Time step
equilibriumAngles = [0, pi];  % Equilibrium positions

% Create figure with subplots
figure;
tiledlayout(2, 2);

% Compare system responses at each equilibrium point
for eqIndex = 1:2
    % Create linear state-space model
    A = A_matrices{eqIndex};
    B = B_matrices{eqIndex};
    linearSys = ss(A, B, C, D);
    
    % Initial condition (slightly perturbed from equilibrium)
    initialAngle = equilibriumAngles(eqIndex) + pi/18;
    initialAngleRate = 0;
    initialState = [initialAngle; initialAngleRate];
    
    % Zero input (free response)
    inputSignal = zeros(size(timeSpan));
    
    % Simulate linear system
    [linearResponse, time, linearStates] = lsim(linearSys, inputSignal, timeSpan, initialState);
    
    % Simulate nonlinear system using Runge-Kutta integration
    nonlinearStates = zeros(length(timeSpan), 2);
    nonlinearStates(1, :) = initialState;
    
    for step = 1:length(timeSpan)-1
        currentTime = timeSpan(step);
        currentState = nonlinearStates(step, :)';
        
        % Fourth-order Runge-Kutta method
        k1 = pendulumDynamics(currentTime, currentState, pendulumMass, cartMass, pendulumLength, gravity, 0);
        k2 = pendulumDynamics(currentTime + timeStep/2, currentState + timeStep/2 * k1, pendulumMass, cartMass, pendulumLength, gravity, 0);
        k3 = pendulumDynamics(currentTime + timeStep/2, currentState + timeStep/2 * k2, pendulumMass, cartMass, pendulumLength, gravity, 0);
        k4 = pendulumDynamics(currentTime + timeStep, currentState + timeStep * k3, pendulumMass, cartMass, pendulumLength, gravity, 0);
        
        % Update state
        nonlinearStates(step+1, :) = (currentState + (timeStep/6) * (k1 + 2*k2 + 2*k3 + k4))';
    end
    
    % Plot nonlinear response
    nexttile;
    plot(time, nonlinearStates(:, 1), 'k', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['Nonlinear Free Response at \theta_{eq} = ', num2str(equilibriumAngles(eqIndex))]);
    grid on;
    
    % Plot linear response
    nexttile;
    plot(time, linearStates(:, 1), 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['Linear Free Response at \theta_{eq} = ', num2str(equilibriumAngles(eqIndex))]);
    grid on;
end

% Nonlinear pendulum dynamics function
function stateDerivative = pendulumDynamics(~, state, m, M, L, g, u)
    theta = state(1);
    thetaDot = state(2);
    
    % Intermediate calculations
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    denominator = M + m * sinTheta^2;
    
    % Nonlinear dynamics equations
    F = (sinTheta / (L * denominator)) * (-m * L * thetaDot^2 * cosTheta + (M + m) * g);
    G = (-cosTheta) / (L * denominator);
    
    thetaDoubleDot = F + G * u;
    
    stateDerivative = [thetaDot; thetaDoubleDot];
end