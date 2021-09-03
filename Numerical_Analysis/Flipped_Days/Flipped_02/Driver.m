f   = @(x) x.^6 - x - 1;
df  = @(x) 6*x.^5 - 1;
ddf = @(x) 30*x.^4;

x0   = 2;
x1   = 1;
root = 1.1347241384015194926;
% root = -0.77808959867860109788
tol  = 1e-15;
maxIts = 100;

% Run each method!
newton(x0, f, df, maxIts, tol, root)

secant(x0, x1, f, maxIts, tol, root)

halley(x0, f, df, ddf, maxIts, tol, root)