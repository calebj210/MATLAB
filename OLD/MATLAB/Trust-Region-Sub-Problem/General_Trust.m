%{
% Solve N-D Elliptical TRS
% Authors: Alan Bouwman, Caleb Jacobs
%}

% Clean workspace
clear
close

% Number of dimensions
n = 200;

% Cost function definitions
% f = @(x) rosenbrock(x);         % Cost function
f = @(x) rastrigin(x);          % Cost function

% Trust region parameters
maxIts  = 500;                   % Max iterations allowed for convergence
epsilon = 10^(-5);               % Stopping tolerance
Delta   = 1;                   % Initial trust region "radius"
scale   = 0.1;                   % Delta changer
x       = ones(n,1);            % Initial minimizer guess

% B = eye(n);

% Begin solving the TRS
for i = 1:maxIts
    % Compute current gradient
    g = finiteGradient(f, x);
    
    % Check stopping criteria
    if (norm(g) < epsilon)
        break;
    end
    
    % Compute current hessian
    % A = H(x);
    A = finiteHessian(f, x);
    
    % Compute Newton step
    p0 = -A \ g;
    
    % Compute elliptical matrix for this step
    B = getEllipticalMatrix(p0, 2);
    
    % Assume we will use Newton step
    useP0 = false;
    
    % See if Newton step is feasible (i.e. ||p0||_B < Delta)
    if BNorm(p0, B) < Delta
        useP0 = true;   % Newton step is feasible, so use it    
    end
    
    % Compute left matrix of M
    M0 = [-B A;A, -g*g'/Delta^2];
          
    % Compute right matrix of M
    M1 = [zeros(n) B;B zeros(n)];
    
    % Find the rightmost eigenvalue/vector of the pencil M_tilde
    [y, lambda] = eigs(M0, -M1, 1, 'largestreal');
    y1 = y(1:n);
    y2 = y(n+1:end);
    
    % Compute boundary step from rightmost eigenvector
    p1 = -sign(g' * y2) * Delta * y1 / BNorm(y1, B);
    
    % Compare Newton step to TR-boundary step to get optimal solution
    if (useP0 && f(x + p0) < f(x + p1))
        x = x + p0;
        fprintf("p0: %d, x err = %.6f, Delta = %.2f\n", i, norm(x - ones(n,1)), Delta)
    else
        x = x + p1;
        fprintf("p1: %d, x err = %.6f, Delta = %.2f\n", i, norm(x - ones(n,1)), Delta)
    end
end

norm(x - ones(n,1))
i

% Elliptical norm
function a = BNorm(p, B)
    a = sqrt(p'*B*p);
end