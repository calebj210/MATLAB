%%
% Problem 2 driver
% Author: Caleb Jacobs
% Date last modified: 09-09-2021
format longE

%% Function parameters
Ti = 20;                % Initial soil temperature [deg C]
Ts = -15;               % Surrounding temperature [deg C]
alpha  = 0.138e-6;      % Thermal conductivety [meters^2/seconds]
t = 5184000;            % Time after cold snap [seconds]

%% Function Definitions
f   = @(x) erf(x / (2 * sqrt(alpha * t))) * (Ti - Ts) + Ts;
df  = @(x) ((Ti - Ts) / sqrt(pi * alpha * t)) * exp(-x.^2 / (4 * alpha * t));

%% Simulation parameters
x0 = 0.01;              % Newton's method initial guess
a = 0;                  % Initial left bound of bisection interval
b = 6;                  % Initial right bound of bisection interval
tol  = 1e-13;           % Stopping tolerance
maxIts = 100;           % Maximum iterations allowed

%% Plot the temperature function
X = linspace(0, 6, 1000);
plot(X, f(X), 'LineWidth', 2)
xlabel('x (meters)')
ylabel('Temperature (deg C)')
title('Temperature of Water Main vs Depth')

%% Run each method!
newton(x0, f, df, maxIts, tol, 0.6769)

bisection(a, b, f, tol, maxIts, 0.6769)