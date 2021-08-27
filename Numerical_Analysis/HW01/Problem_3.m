% Homework 1 problem 3
% Author: Caleb Jacobs
% Date last modified 26-08-2021

% Function defintions
f1 = @(x, delta) cos(x + delta) - cos(x);
f2 = @(x, delta) -2 * sin((2*x + delta) / 2) .* sin(delta / 2);
f3 = @(x, delta, xi) -delta * sin(x) - (delta .^ 2) / 2 * cos(xi);

% Variable defintions
X = 0;
Deltas = logspace(-16, 0);

% Plot differences
loglog(Deltas, abs(f1(X, Deltas) - f2(X, Deltas)), 'LineWidth', 1.5)

% loglog(Deltas, abs(f1(X, Deltas) - f3(X, Deltas, X)), 'LineWidth', 1.5)