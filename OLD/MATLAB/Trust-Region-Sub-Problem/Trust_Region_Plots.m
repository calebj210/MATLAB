%{
% Solve N-D Elliptical TRS
% Authors: Alan Bouwman, Caleb Jacobs
%}

% Clean workspace
clear
close

% Number of dimensions
n = 2;

% Cost function definitions
% Cost function
f = @(x) rosenbrock(x);          % Cost function
% f = @(x) rastrigin(x);         % Cost function

% Trust region parameters
maxIts  = 200;                   % Max iterations allowed for convergence
epsilon = 10^(-5);               % Stopping tolerance
Delta   = 0.5;                  % Initial trust region "radius"
scale   = 0.1;                   % Delta changer
x       = [-5;5];                % Initial minimizer guess

% Elliptical trust region matrix
va = [1;1];                      % Major axis vector
vb = [1;-1];                     % Minor axis vector
a = 1;                           % Radius of major axis
b = 1;                           % Radius of minor axis

% B = getEllipticalMatrix(a, va, b, vb);

% Compute Rosenbrock over designated area
figure();
[X, Y] = meshgrid(-5:0.1:2, -2:0.1:7);
F = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        F(i, j) = f([X(i,j), Y(i,j)]);
    end
end

% Create contour plot of Rosenbrock
contour(X, Y, F, 200)
hold on
colorbar
title('Elliptical Rosenbrock Optimization Path')
pbaspect([7 9 1])

% Begin solving the TRS
for i = 0:maxIts
    % Compute current gradient
    g = finiteGradient(f, x);
    
    % Check stopping criteria
    if (norm(g) < epsilon)
        break;
    end
    
    % Compute current hessian
    A = finiteHessian(f, x);
    
    % Compute Newton step
    p0 = -A \ g;
    
    % Compute elliptical matrix
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
        fprintf("p0: %d, x = (%.1f,%.1f), Delta = %.2f\n", i, x(1), x(2), Delta)
    else
        x = x + p1;
        fprintf("p1: %d, x = (%.1f,%.1f), Delta = %.2f\n", i, x(1), x(2), Delta)
    end
    
    % Add the new iteration to the plot
    scatter(x(1), x(2), 200, '.')
    
    % Compute trust region boundaries
    t = 0 : 0.01 : 2 * pi;
    r = zeros(length(t), 2);
    [V,D] = eig(B);
    va = V(:,1) / norm(V(:,1));
    vb = V(:,2) / norm(V(:,2));
    a = D(1,1);
    b = D(2,2);
    for j = 1 : length(t)
        r(j,:) = x + Delta * (a * va * cos(t(j)) + b * vb * sin(t(j)));
    end
    
    % Display trust region boundary
    plot(r(:,1), r(:,2))
end

B
x
i

% Elliptical norm
function a = BNorm(p, B)
    a = sqrt(p'*B*p);
end

% % Generate elliptical region matrix
% function B = getEllipticalMatrix(a, va, b, vb)
%     va = va / norm(va);     % Normalize major axis vector
%     vb = vb / norm(vb);     % Normalize minor axis vector
%     
%     P = [va vb];
%     D = diag([1/a,1/b]).^2;
%     
%     B = P * D * P^(-1);
% end