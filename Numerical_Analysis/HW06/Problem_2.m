%%
% Homework 6, problem 2 code
% Linear system solving using SD and CG
%
% Author: Caleb Jacobs
% Date last modified: 6-10-2021

close all
clear
format longE

%% Settings
tau    = [0.01; 0.05; 0.1; 0.2];  % Cutoff values
N      = 500;                     % Size of matrices
seed   = 210;                     % Random number seed
tol    = 1e-10;                   % Error Toloerance
maxIts = 50;                      % Maximum allowed iterations
x0     = zeros(N, 1);             % Initial solution guess


%% Linear system setup
% Construct test matrices
A = zeros(N, N, length(tau));
for i = 1 : length(tau)
    A(:,:,i) = genMat(N, tau(i), seed);
end

% Construct b
b = rand(N, 1);

% True solutions for error calculations
xTrue = zeros(N, length(tau));
for i = 1 : length(tau)
    xTrue(:, i) = A(:, :, i) \ b;
end

%% Driver
% Run steepest descent
figure()
x1 = zeros(N, length(tau));
for i = 1 : length(tau)
    [x1(:, i), r] = sd(A(:,:,i), b, x0, tol, maxIts);
    
%     norm(A(:,:,i) * x1(:,i) - b);
    semilogy(r, 'LineWidth', 2)
    hold on
end
title('Steepest Descent Norm(Residual) vs Iterations')
xlabel('Iterations')
ylabel('||r_k||')
legend('\tau = 0.01', '\tau = 0.05', '\tau = 0.1', '\tau = 0.2')

% Run conjugate gradient
figure()
x2 = zeros(N, length(tau));
for i = 1 : length(tau)
    [x2(:, i), r] = cg(A(:,:,i), b, x0, tol, maxIts);
    
%     norm(A(:,:,i) * x2(:,i) - b);
    semilogy(r, 'LineWidth', 2)
    hold on
end
title('Conjugate Gradient Norm(Residual) vs Iterations')
xlabel('Iterations k')
ylabel('||r_k||')
legend('\tau = 0.01', '\tau = 0.05', '\tau = 0.1', '\tau = 0.2')

%% Compute error bounds
% Get largest and smallest eigenvalues of each matrix
lambda = zeros(length(tau), 2);
for i = 1 : length(tau)
    lambda(i, 1) = eigs(A(:,:,i), 1, 'smallestabs');   % Smallest modulus eigenvalue
    lambda(i, 2) = eigs(A(:,:,i), 1, 'largestabs');    % Largest modulus eigenvalue
end

% Get condition number of each matrix
conds = zeros(length(tau), 1);
for i = 1 : length(tau)
    conds(i) = cond(A(:,:,i));
end

% SD error coefficients
sdErrCoef = zeros(length(tau), 1);
for i = 1 : length(tau)
    sdErrCoef(i) = sqrt((lambda(i, 2) - lambda(i, 1)) ./ (lambda(i, 2) + lambda(i, 1)));
end
sdErrCoef

% CG error bounds
c = (sqrt(conds) - 1) ./ (sqrt(conds) + 1)

%% Generate random matrix with specified tau
function A = genMat(n, tau, s)
    rng(s);                             % Reset random # generator
    tmp = 2 * rand(n-1, 1) - 1;         % Random vector
    
    % Strip random vector
    for i = 1 : n-1
        if abs(tmp(i)) > tau
            tmp(i) = 0;
        end
    end
    
    A = diag(tmp, -1) + diag(ones(n, 1)) + diag(tmp, 1); 
end

%% Steepest Descent solver
function [x, r] = sd(A, b, x0, tol, maxIts)
    rk = 1 - A * x0;            % Initialize residual
    
    r = zeros(maxIts + 1, 1);   % Initial residual norm storage
    r(1) = norm(rk);            % Store first residual
    
    for i = 1 : maxIts
        % Check convergence criteria
        if r(i) < tol
            break;
        end
        
        % Compute next iterate
        x0 = x0 + (rk' * rk) / (rk' * A * rk) * rk;
        
        rk = b - A * x0;        % Compute new residual
        r(i + 1) = norm(rk);    % Store new residual
    end
    
    x = x0;                     % Return solution
end

%% Conjugate Gradient solver
function [x, r] = cg(A, b, x0, tol, maxIts)
    rk  = b - A * x0;                   % Initial residual
    p   = rk;                           % Initial search direction
    rri = rk' * rk;                     % ||r_k||^2
    
    r    = zeros(maxIts + 1, 1);        % Full residual storage
    r(1) = sqrt(rri);                   % Store initial residual
    
    for i = 1:maxIts
        Ap  = A * p;                    % A*p
        a   = rri / (p' * Ap);          % search length
        x0  = x0 + a * p;               % Get next iterate
        rk  = rk - a * Ap;              % Compute new residual
        rrf = rk' * rk;                 % Compute new ||r_k||^2
        
        r(i + 1) = sqrt(rrf);           % Store new ||r_k||
        
        % Check convergence criteria
        if sqrt(rrf) < tol
            break;
        end
        
        p = rk + (rrf / rri) * p;       % Compute next search direction
        rri = rrf;                      % Save ||r_k||^2
    end
    
    x = x0;                             % Return solution
end