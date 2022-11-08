%%
% APPM 5600 Iterative Solvers
% Author: Caleb Jacobs
% Date Last Modified: 01-10-2021

%% Problem parameters
A = [4,-1,0,-1,0,0;  ...
     -1,4,-1,0,-1,0; ...
     0,-1,4,-1,0,-1; ...
     -1,0,-1,4,-1,0; ...
     0,-1,0,-1,4,-1; ...
     0,0,-1,0,-1,4];
b = [2;1;2;2;1;2];

%% Settings
format long
x0 = [1;1;1;1;1;1];
maxIts = 100;
tol = 1e-7;
omega = 1.6735;

%% Driver for iterative methods
GJ(x0, A, b, tol, maxIts)
GS(x0, A, b, tol, maxIts)
SOR(x0, A, b, omega, tol, maxIts)

%% Question 3 Spectral Radii
D = diag(diag(A));  % Diagonal entries of A
L = tril(A,-1);      % Lower triangular entries of A
U = triu(A,1);      % Upper triangular entries of A
jacSpec  = abs(eigs(-D\(L + U), 1, 'largestabs'))
seidSpec = abs(eigs(-(L + D) \ U, 1, 'largestabs'))
SORSpec  = abs(eigs(-(D + omega * L) \ (omega * U + (omega - 1) * D), 1, 'largestabs'))

%% Question 4 Spectral Radii
omega = 0.8 : 0.1 : 1.3;
A = @(w) [1 - w, w / 2; w .* (1 - w) / 2, (w - 2).^2 / 4];
spec = zeros(length(omega), 1);
for i = 1:length(omega)
    spec(i) = eigs(A(omega(i)), 1, 'largestabs');
end
[sig, idx] = min(abs(spec));
sig                         % Minimum spectral radius
omega(idx)                  % Omega spectrum minimizer

%% Gauss-Jacobi iterative solver
function x = GJ(x0, A, b, tol, maxIts)
    fprintf('Gauss-Jacobi Iteration\n')
    k = 0;
    xf = x0;
    n = length(x0);
    
    
    while 1
        xi = xf;
        
        for i = 1:n
            tmp = 0;
            
            for j = 1:n
                if j ~= i
                    tmp = tmp + A(i,j) * xi(j);
                end
            end
            
            xf(i) = (b(i) - tmp) / A(i,i);
        end
        
        if norm(xf - xi,inf) < tol
            x = xf;
            k
            break;
        end
        if k >= maxIts
            x = NaN;
            break;
        end
        k = k + 1;
    end
end

%% Gauss-Seidel iterative solver
function x = GS(x0, A, b, tol, maxIts)
    fprintf('Gauss-Siedel Iteration\n')
    k = 0;
    xf = x0;
    n = length(x0);
    
    
    while 1
        xi = xf;
        
        for i = 1:n
            tmp = 0;
            
            for j = 1 : i-1
                tmp = tmp + A(i,j) * xf(j);
            end
            
            for j = i+1 : n
                tmp = tmp + A(i,j) * xf(j);
            end
            
            xf(i) = (b(i) - tmp) / A(i,i);
        end
        
        if norm(xf - xi, inf) < tol
            x = xf;
            k
            break;
        end
        if k >= maxIts
            x = NaN;
            break;
        end
        k = k + 1;
    end
end

%% SOR Iterative Solver
function x = SOR(x0, A, b, omega, tol, maxIts)
    fprintf('SOR Iteration\n')
    k = 0;
    xf = x0;
    n = length(x0);
    
    while 1
        xi = xf;
        
        for i = 1:n
            tmp = 0;
            
            for j = 1:n
                if j ~= i
                    tmp = tmp + A(i,j) * xf(j);
                end
            end
            
            xf(i) = (1 - omega) * xf(i) + omega * (b(i) - tmp) / A(i,i);
        end
        
        if norm(xf - xi, inf) < tol
            x = xf;
            k
            break;
        end
        if k >= maxIts
            x = NaN;
            break;
        end
        k = k + 1;
    end
end