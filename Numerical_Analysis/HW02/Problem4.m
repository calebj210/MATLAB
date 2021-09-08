%% Homework 2 Newton's method and Modified Newton
% Author: Caleb Jacobs
% Date last modified: 08-09-2021

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
newton(x0, f, df, maxIts, 1e-16, root)
modNewton(x0, f, df, p, maxIts, tol, root)

%% Function defintions
% Newton's method
function ier = newton(x0, f, df, maxIts, tol, root)
    fprintf("Newton's Method\n")
    for i = 1:maxIts
        denom = df(x0);
        
        if denom == 0
            if f(x0) == 0
                ier = 0;
                fprintf('Root = %.16f\n\n', x1)
            else
                ier = 2;
            end
            return
        end
        
        x1 = x0 - f(x0) / denom;
        
        if abs(x1 - x0) <= tol
            ier = 0;
            fprintf('Root = %.16f\n\n', x1);
            return
        end
        
        x0 = x1;
        
        fprintf('%3d: %.16e\n', i, abs(x1 - root) / root)
    end
    
    ier = 1;
    return
end

% Modfied Newton's method
function ier = modNewton(x0, f, df, p, maxIts, tol, root)
    fprintf("Modified Newton's Method\n")
    for i = 1:maxIts
        denom = df(x0);
        
        if denom == 0
            if f(x1) == 0
                ier = 0;
                fprintf('Root = %.16f\n\n', x1)
            else
                ier = 2;
            end
            return
        end
        
        x1 = x0 - p * f(x0) / denom;
        
        if abs(x1 - x0) <= tol
            ier = 0;
            fprintf('Root = %.16f\n\n', x1);
            return
        end
        
        x0 = x1;
        
        fprintf('%3d: %.16e\n', i, abs(x1 - root) / root)
    end
    
    ier = 1;
    return
end