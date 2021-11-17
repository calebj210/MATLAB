%%
% Composite trapezoidal and Simpson's rule convergence
%
% Author: Caleb Jacobs
% Date last modified: 17-Nov-2021

%% Settings
format long

%% Parameters
f = @(x) -4 * x .* log(x);      % Function to integrate
a = 0;                          % Left endpoint
b = 1;                          % Right endpoint 
N = 2 .^ [1 : 9];               % Create interval array

%% Compute integrals
errs = zeros(length(N), 3);
for i = 1 : length(N)
    n = N(i);
    EM = abs(mid(f, a, b, n) - 1);
    ET = abs(trapz(f, a, b, n) - 1);
    ES = abs(simps(f, a, b, n) - 1);
    
    errs(i, :) = [EM, ET, ES];
end

figure(1)
loglog((b - a) ./ N, errs, 'LineWidth', 2)
legend('Midpoint', 'Trapezoidal', "Simpson's")
title('Convergence of Different Quadrature')
xlabel('Step size (h)')
ylabel('Absolute Error')

%% Function defintions
% Composite trapezoidal rule
function val = trapz(f, a, b, n)
    xi = linspace(a, b, n + 1);         % Compute evaluation points
    h = xi(2) - a;                      % Compute x spacing
    
    if ~isnan(f(a))
        val = h * (f(a) + f(b)) / 2;        % Add endpoint contribution
    else
        val = h * f(b) / 2;
    end
    val = val + h * sum(f(xi(2 : n)));  % Add interior contribtuion
end

% Composite Simposon's rule
function val = simps(f, a, b, n)
    % Round n up to nearest even number
    if mod(n, 2) == 0
        N = n;
    else
        N = n + 1;
    end
    
    xi = linspace(a, b, N + 1);             % Compute evaluation points
    h = xi(2) - a;                          % Compute x spacing
    
    if isnan(f(a))
        val = f(b);                         % Add endpoint contribution
    else
        val = f(a) + f(b);
    end
    
    val = val + 4 * sum(f(xi(2 : 2 : N)));  % Add odd node contribution
    
    val = val + 2 * sum(f(xi(3 : 2 : N)));  % Add even node contribution
    
    val = h * val / 3;                      % Scale integral accordingly
end

% Composite midpoint rule
function val = mid(f, a, b, n)
    xTmp = linspace(a, b, n + 1);               % Compute standard points
    xi   = (xTmp(2 : n + 1) + xTmp(1 : n)) / 2; % Get midpionts
    h    = xTmp(2) - a;                         % Compute step size
    
    val = h * sum(f(xi));                       % Compute integral
end