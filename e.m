clc; clear; close all;

% System parameters
pendulumMass = 0.1;     % Mass of the pendulum (kg)
cartMass = 1;           % Mass of the cart (kg)
pendulumLength = 0.5;   % Length of the pendulum (m)
gravity = 9.81;         % Acceleration due to gravity (m/s²)

% Linearized state-space matrices for two equilibrium points
A = {
    [0, 1; (cartMass + pendulumMass) * gravity / (pendulumLength * cartMass), 0],  % Around θ = 0
    [0, 1; -(cartMass + pendulumMass) * gravity / (pendulumLength * cartMass), 0]  % Around θ = π
};
B = {
    [0; -1/(pendulumLength * cartMass)],  % Input matrix for θ = 0
    [0; 1/(pendulumLength * cartMass)]    % Input matrix for θ = π
};

% Output and feedthrough matrices
C = [1, 0];  % Measure angular position
D = 0;       % No direct feedthrough

% Time vector for simulation
timeSpan = linspace(0, 20, 500);
initialConditions = [0, 0; pi, 0];  % [θ, dθ/dt] for each equilibrium

% Create figure with subplots
figure;
tiledlayout(2, 2);

for eqPoint = 1:2
    % Linear system response
    linearSys = ss(A{eqPoint}, B{eqPoint}, C, D);
    [linearTheta, time, linearStates] = lsim(linearSys, ones(size(timeSpan)), timeSpan, initialConditions(eqPoint, :));
    
    % Nonlinear system response using ODE solver
    [timeNonlin, statesNonlin] = ode45(@(t, y) pendulumDynamics(t, y, pendulumMass, cartMass, pendulumLength, gravity, 1), timeSpan, initialConditions(eqPoint, :));
    
    % Plot nonlinear response
    nexttile;
    plot(timeNonlin, statesNonlin(:, 1), 'k', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['Nonlinear Step Response at \theta_{eq} = ', num2str(initialConditions(eqPoint, 1))]);
    grid on;
    
    % Plot linear response
    nexttile;
    plot(time, linearTheta, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['Linear Step Response at \theta_{eq} = ', num2str(initialConditions(eqPoint, 1))]);
    grid on;
end

% Nonlinear pendulum dynamics
function dydt = pendulumDynamics(~, state, m, M, L, g, u)
    theta = state(1);
    thetaDot = state(2);
    
    % Nonlinear dynamics equations
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    denominator = M + m * sinTheta^2;
    
    thetaDDot = (sinTheta/ (L * denominator)) * (-m * L * thetaDot^2 * cosTheta + (M + m) * g) ...
                - (cosTheta / (L * denominator)) * u;
    
    dydt = [thetaDot; thetaDDot];
end