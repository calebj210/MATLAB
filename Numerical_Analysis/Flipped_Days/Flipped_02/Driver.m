f   = @(x) (x-1).*(x-2).*(x-3);
df  = @(x) (x-1).*(x-2) + (x-2).*(x-3) + (x-2).*(x-3)
ddf = @(x) 6*x;

x0   = 0.0000000001;
x1   = 1;
root = 0;
% root = -0.77808959867860109788;
tol  = 1e-15;
maxIts = 100;

% Run each method!
newton(x0, f, df, maxIts, tol, root)

% secant(x0, x1, f, maxIts, tol, root)

% halley(x0, f, df, ddf, maxIts, tol, root)