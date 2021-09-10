Ti = 20;
Ts = -15;
alpha  = 0.138e-6;
t  = 5184000;

f   = @(x) erf(x / (2 * sqrt(alpha * t))) * (Ti - Ts) + Ts;
df  = @(x) ((Ti - Ts) / sqrt(pi * alpha * t)) * exp(-x.^2 / (4 * alpha * t));
ddf = @(x) 6*x;

X = linspace(0, 6, 1000);
plot(X, f(X), 'LineWidth', 2)
xlabel('x (meters)')
ylabel('Temperature (deg C)')
title('Temperature of Water Main vs Depth')

x0   = 0.01;
x1   = 1;
a = 0;
b = 6;
root = 0;
% root = -0.77808959867860109788;
tol  = 1e-15;
maxIts = 100;

% Run each method!
newton(x0, f, df, maxIts, tol, root)

bisection(a, b, f, tol, maxIts)

% secant(x0, x1, f, maxIts, tol, root)

% halley(x0, f, df, ddf, maxIts, tol, root)