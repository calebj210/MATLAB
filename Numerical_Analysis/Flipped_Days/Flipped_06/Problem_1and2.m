%%
% Trig interpolation routine with FFT
% Author: Caleb Jacobs
% Date last modified: 21-10-2021

%% Intitial setup
format longe
clear

%% Function and data setup
f = @(x) exp(sin(x));           % Function to interpolate

a = 0;                          % Left boundary
b = 2*pi;                       % Right boundary
N = 10;                         % Number of nodes

x = linspace(a, b, N);          % Node set
y = f(x);                       % Data set

t = linspace(a, b, 1000)';      % Evaluation set

%% Driver
p = interp(t, y);

figure(1)
plot(t, p, t, f(t))

%% FFT interpolant
function F = interp(t, y)
    n = length(y);
    m = length(t);
    c = fft(y) / n;
    k = [0 : n - 1]';
    
    F = zeros(m, 1);
    for i = 1 : m
        F(i) = c * exp(1i * t(i) * k);
    end
end