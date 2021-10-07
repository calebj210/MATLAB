%%
% Homework 6, problem 2 code
% Author: Caleb Jacobs
% Date last modified: 6-10-2021
close all
clear

%% Settings
tau    = [0.01; 0.05; 0.1; 0.2];  % Cutoff values
N      = 500;                     % Size of matrices
seed   = 210;                     % Random number seed
tol    = 1e-10;                   % Error Toloerance
maxIts = 1000;                    % Maximum allowed iterations
x0     = zeros(N, 1);             % Initial solution guess


%% Linear system setup
% Construct test matrices
A = zeros(N, N, length(tau));
for i = 1 : length(tau)
    A(:,:,i) = genMat(N, tau(i), seed);
end

% Construct b
b = rand(N, 1);


%% Driver
x1 = zeros(N, length(tau));
for i = 1 : length(tau)
    [x1(:, i), r] = sd(A(:,:,i), b, x0, tol, maxIts);
    
    norm(A(:,:,i) * x1(:,i) - b)
    semilogy(r)
    hold on
end

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
function x = cg(A, b, x0, tol, maxIts)
    rk = b - A * x0;            % Initialize residual
    pk = rk;                    % Initial conjugate vector
    
    r = zeros(maxIts + 1, 1);   % Initial residual norm storage
    r(1) = norm(rk);            % Store first residual
    
    for i = 1 : maxIts
        ak = (rk' * rk) / (pk' * A * pk);   % Compute distance
        
        x0 = x0 + ak * pk;          % Next iterate
        rkNew = rk - ak * A * x0;   % Compute new residual
        
        % Check convergence criteria
        if rkNew < tol
            break;
        end
        
        r(i + 1) = norm(rk);    % Store new residual
    end
    
    x = x0;                     % Return solution
end