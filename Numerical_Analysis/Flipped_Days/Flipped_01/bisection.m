% Bisection routine for computing roots of functions
% Author: Caleb Jacobs
% Date last modified: 02-09-2021

format longE

% Function defintions
f1 = @(x) (x - 5).^9;
f2 = @(x) -1953125 + 3515625*x - 2812500*x.^2 + 1312500*x.^3 - ...
          393750*x.^4 + 78750*x.^5 - 10500*x.^6 + ...
          900*x.^7 - 45*x.^8 + x.^9;

% Set parameters
a   = 4.8;
b   = 5.31;
tol = 1e-4;
maxIts = 20;

% Find root with bisection
bisect(a, b, f1, tol, maxIts);
bisect(a, b, f2, tol, maxIts);

% Bisection method definition
function c = bisect(a, b, f, tol,  maxIts)
    % Check if root is guaranteed in the interval
    if f(a) * f(b) > 0
        fprintf('Root not guarnteed, exiting!\n\n')
        return
    end
    
    for i = 1:maxIts
        c = (b + a) / 2;                    % Midpoint of interval
        
        if f(c) == 0 || (b - c) < tol       % Check convergence criteria
            fprintf('Its used = %d\n',i)    % Display number of iterations
            fprintf('x = %.16f\n\n', c)     % Display root information
            return
        end
        
        if f(c) * f(a) >= 0                 % Update searching interval
            a = c;
        else
            b = c;
        end
    end
    
    fprintf('Convergence could not happen to desired tolerance\n\n')
end