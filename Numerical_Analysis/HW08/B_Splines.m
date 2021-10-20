%%
% B-Spline interpolation routine for Problem 3 of HW08
% Author: Caleb Jacobs
% Date last modified: 20-10-2021

%% Clean workspace
clear

%% Settings
format longE
a = 0;          % Left side of interval
b = 2*pi;       % Right side of interval
n = 20;          % Number of intervals
m = 25;         % Number of data points

%% Initialize data
h  = (b - a) / m;               % Step size
x  = linspace(a, b, n + 1);     % Interval boundaries
t  = linspace(a, b, m);         % Nodes

A = constructA(t, x);           % Find A matrix for computing weights

%% Interpololate sin(x)
f = @(x) sin(x);                % Function definition
bb = f(t)' + 0.05*randn(m, 1);  % Data function values
c = A \ bb;                     % Interpolation weights

figure(1)
T = linspace(a, b, 1000);
plot(T, S(T, x, c), T, f(T), 'LineWidth', 1.5)
title('sin(x) Interpolant with Noisy Data')
xlabel('x')
ylabel('f(x)')

%% Interpolate exp(x)
f = @(x) exp(x);                % Function definition
bb = f(t)' + 6*randn(m, 1);     % Data function values
c = A \ bb;                     % Interpolation weights

figure(2)
T = linspace(a, b, 1000);
plot(T, S(T, x, c), T, f(T), 'LineWidth', 1.5)
title('e^x Interpolant with Noisy Data')
xlabel('x')
ylabel('f(x)')

%% Interpolate sin(x)/x
f = @(x) sin(4*x) ./ (4*x + 1); % Function definition
bb = f(t)' + 0.02*randn(m, 1);  % Data function values
c = A \ bb;                     % Interpolation weights

figure(3)
T = linspace(a, b, 1000);
plot(T, S(T, x, c), T, f(T), 'LineWidth', 1.5)
title('sin(4x)/(4x+1) Interpolant with Noisy Data')
xlabel('x')
ylabel('f(x)')

%% Plot matrix structure
figure(4)
spy(A)
title('Spy plot of A matrix')
xlabel('Columns')
ylabel('Rows')

%% Routine definitions
% Basis spline
function val = B(t, x, i)
    h = x(2) - x(1);
    n = length(x);
    
    if i - 2 > 0 && t < x(i - 2)
        val = 0;
    elseif i - 2 > 0 && x(i - 2) <= t && t <= x(i - 1)
        val = (t - x(i - 2)).^3 / h^3;
    elseif i - 1 > 0 && x(i - 1) <= t && t <= x(i)
        xi  = x(i - 1);
        val = 1 + 3*(t - xi) / h + 3*(t - xi).^2 / h^2 - 3*(t - xi)^3 / h^3;
    elseif i + 1 <= n && x(i) <= t && t <= x(i + 1)
        xi  = x(i + 1);
        val = 1 + 3*(xi - t) / h + 3*(xi - t).^2 / h^2 - 3*(xi - t)^3 / h^3;
    elseif i + 2 <= n && x(i + 1) <= t && t <= x(i + 2)
        val = (x(i + 2) - t).^3 / h^3;
    else
        val = 0;
    end
end

% Get interval indices of data
function idx = getIdx(t, x)
    m = length(t);          % Number of nodes
    n = length(x);          % Number of intervals + 1
    
    idx = ones(m, 1);
    for i = 1 : m
        for j = 2 : n
            if t(i) < x(j)
                break
            end
        end
        idx(i) = j - 1;
    end
end

% Constuct A matrix
function A = constructA(t, x)
    m = length(t);                  % Number of nodes
    n = length(x);                  % Number of intervals + 1
    h = x(2) - x(1);                % Data spacing
    
    idx = getIdx(t, x) + 2;         % Interval indices of data
   
    x = [(x(1) - h*[2 1]) x (x(n) + h*[1 2])];  % Append extra intervals
    
    A = spalloc(m, n + 2, 4 * m);   % Initialize sparse A
    for i = 1 : m
        for j = -1 : 2
            A(i, idx(i) + j) = B(t(i), x, idx(i) + j);
        end
    end
end

% Evaluate Interpolant
function f = S(t, x, c)
    m = length(t);
    h = x(2) - x(1);
    n = length(x);
    x = [(x(1) - h*[2 1]) x (x(n) + h*[1 2])];  % Append extra intervals
    n = length(x);
    
    f = zeros(m, 1);
    for j = 1 : m
        for i = 1 : n - 1
            f(j) = f(j) + c(i) * B(t(j), x, i);
        end 
    end
end