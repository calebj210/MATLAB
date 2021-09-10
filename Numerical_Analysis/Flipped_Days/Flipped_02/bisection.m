% Bisection method definition
function c = bisection(a, b, f, tol,  maxIts)
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