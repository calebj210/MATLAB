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
n = 10;                         % Number of nodes
N = 1000;                       % Number of evaluation nodes

x = linspace(a, b, n);          % Node set
y = f(x);                       % Data set

t = linspace(a, b, N)';         % Evaluation set

%% Driver
p1 = interp(t, y)';
p2 = interpft(y, N);
% p = abs(ifft(fft(y) * N / n, N));

figure(1)
plot(t, p1, t, p2, t, f(t))

%% FFT interpolant
function F = interp(t, y)
    n = length(y);
    m = length(t);
    c = fftshift(fft(y) / n);
    k = [0 : n - 1]';
    
    F = zeros(m, 1);
    for i = 1 : m
        F(i) = abs(c * exp(1i * t(i) * k));
    end
end