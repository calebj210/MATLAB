%%
% B-Spline interpolation routine for Problem 3 of HW08
% Author: Caleb Jacobs
% Date last modified: 18-10-2021

%% Clean workspace
clear

%% Settings
format longE
a = 0;          % Left side of interval
b = 1;          % Right side of interval
n = 10;         % Number of intervals
m = 15;         % Number of data points

%% Function definition
f = @(x) sin(x);

%% Initialize data
h  = (b - a) / m;               % Step size
x  = linspace(a, b, n + 1);     % Interval boundaries
t  = linspace(a, b, m);         % Nodes
fi = f(t);                      % Date nodes

%% Routine definitions
% Basis spline
function val = B(t, x, idx, i)
    
end

% Get interval indices of data
function idx = getIdx(t, x)
    m = length(t);          % Number of nodes
    n = length(xi);         % Number of intervals + 1
    
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
    m = length(t);          % Number of nodes
    n = length(xi);         % Number of intervals + 1
    h = t(2) - t(1);        % Data spacing
    
    idx = getIdx(t, x);     % Interval indices of data
    
    A = sparse(m, n + 2);   % Initialize sparse A
    for i = 1 : m
        for j = -1 : 2
            A(i, idx(i) + j) = B(t(i), x, idx, idx(i) + j);
        end
    end
end