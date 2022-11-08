%% Homework 2 Newton's method and Modified Newton
% Author: Caleb Jacobs
% Date last modified: 08-09-2021
format longE
close  all

%% Function definitions
p  = 5;
f  = @(x) (x - 1).^p * exp(x);
df = @(x) p*(x - 1).^(p - 1) * exp(x) + (x - 1).^5 * exp(x);

%% Iteration settings
x0     = 0;         % Initial root guess
tol    = 1e-9;      % Accuracy goal
maxIts = 1000;      % Maximum allowed iterations
root   = 1;         % True root to our function

%% Run root finding functions
dat1 = newton(x0, f, df, maxIts, 1e-16, root);
dat2 = modNewton(x0, f, df, p, maxIts, tol, root);

%% Plot convergence data
figure()
semilogy(0:(length(dat1) - 1), dat1(:,2), 'LineWidth', 2)
title("Newton's Method: Error vs Iterations")
xlabel('Iteration #')
ylabel('Absolute error')

figure()
semilogy(0:(length(dat2) - 1), dat2(:,2), 'LineWidth', 2)
title("Modified Newton's Method: Error vs Iterations")
xlabel('Iteration #')
ylabel('Absolute error')

%% Function defintions
% Newton's method
% x0: initial root guess
% f:  function to find root of
% df: derivative of function 
% maxIts: maximum iterations for convergence
% tol: convergence tolerance
% root: true root for error computation
% return: iteration data
function data = newton(x0, f, df, maxIts, tol, root)
    data = [x0 abs(x0 - root)];     % Initialize data set
    
    % Begin running iterations
    for i = 1:maxIts
        denom = df(x0);             % Compute derivative of f
        
        if denom == 0               % Check if derivative is flat
            fprintf('ERROR: Derivative is flat\n');
            return
        end
        
        x1 = x0 - f(x0) / denom;    % Compute next iteration
        
        data(i + 1, :) = [x1 abs(x1 - root)];   % Add data to array
        
        if abs(x1 - x0) <= tol      % Check stopping criteria
            return
        end
        
        x0 = x1;                    % Save current iteration
    end
    
    fprintf('ERROR: Method did not converge\n')
    return
end

% Modfied Newton's method
% x0: initial root guess
% f:  function to find root of
% df: derivative of function
% p: Multiplicity of root
% maxIts: maximum iterations for convergence
% tol: convergence tolerance
% root: true root for error computation
% return: iteration data
function data = modNewton(x0, f, df, p, maxIts, tol, root)
    data = [x0 abs(x0 - root)];     % Initialize data set
    
    % Begin running iterations
    for i = 1:maxIts
        denom = df(x0);             % Compute derivative of f
        
        if denom == 0               % Check if derivative is flat
            fprintf('ERROR: Derivative is flat\n');
            return
        end
        
        x1 = x0 - p * f(x0) / denom;    % Compute next iteration
        
        data(i + 1, :) = [x1 abs(x1 - root)];   % Add data to array
        
        if abs(x1 - x0) <= tol      % Check stopping criteria
            return
        end
        
        x0 = x1;                    % Save current iteration
    end
    
    fprintf('ERROR: Method did not converge\n')
    return
end