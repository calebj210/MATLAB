% Bisection routine for computing roots of functions
% Author: Caleb Jacobs
% Date last modified: 26-08-2021

format longE

% Functions to find the roots of
f1 = @(x) (x - 1) * (x - 3) * (x - 5);
f2 = @(x) (x - 1)^2 * (x - 3);
f3 = @(x) sin(x);

% Compute the roots of the functions above
bisect(0, 2.4, f1, 10^(-5), 100)        % f1 root finding

bisect(0, 2, f2, 10^(-5), 100)          % f2 root finding

bisect(0, 0.1, f3, 10^(-5), 100)        % f3 root finding
bisect(0.5, 3*pi/4, f3, 10^(-5), 100)   % f3 root finding


% Bisection algorithm definition
function c = bisect(a, b, f, tol,  maxIts)
    i = 1;
    
    while i <= maxIts
        c = (b + a) / 2;                    % Midpoint of interval
        
        if f(c) == 0 || (b - a) / 2 < tol   % Check convergence criteria
            i                               % Display number of iterations
            break;
        end
        
        i = i + 1;
        
        if sign(f(c)) == sign(f(a))         % Update searching interval
            a = c;
        else
            b = c;
        end
    end
end