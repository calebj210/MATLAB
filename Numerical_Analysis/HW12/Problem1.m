%%
% Composite trapezoidal and Simpson's rule timings
%
% Author: Caleb Jacobs
% Date last modified: 17-Nov-2021

%% Settings
format long

%% Parameters
f = @(x) 1 ./ (1 + x.^2);       % Function to integrate
a = -5;                         % Left endpoint
b =  5;                         % Right endpoint
n1 = 108;                       % Number of intervals for Simpson's
n2 = 1291;                      % Number of intervals for Trapezoidal

%% Compute integrals and timings
fprintf('Trapezoidal runtime\n')
tic 
    T = trapz(f, a, b, n2);
toc

fprintf("\nSimpson's runtime\n")
tic
    S = simps(f, a, b, n1);
toc

fprintf('\nMidpoint runtime\n')
tic
    M = mid(f, a, b, n2);
toc

fprintf('\nQuad runtime with error of 10^-4\n')
tic
    [Q1, nQ1] = quad(f, a, b, 1e-4);
toc

fprintf('\nQuad runtime with error of 10^-6\n')
tic
    [Q2, nQ2] = quad(f, a, b, 1e-6);
toc

%% Display quadrature results
fprintf('\nTrapezoidal\tf evals = %d\t%.10f\n', n2 + 1, T)
fprintf("Simpson's\tf evals = %d\t%.10f\n", n1 + 1, S)
fprintf('Midpoint\tf evals = %d\t%.10f\n', n2, M)
fprintf('Quad 10^-4\tf evals = %d\t%.10f\n', nQ1, Q1)
fprintf('Quad 10^-6\tf evals = %d\t%.10f\n', nQ2, Q2)

%% Function defintions
% Composite trapezoidal rule
function val = trapz(f, a, b, n)
    xi = linspace(a, b, n + 1);         % Compute evaluation points
    h = xi(2) - a;                      % Compute x spacing
    
    val = h * (f(a) + f(b)) / 2;        % Add endpoint contribution
    
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
    
    val = f(a) + f(b);                      % Add endpoint contribution
    
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