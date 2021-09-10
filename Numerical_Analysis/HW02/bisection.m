%%
% Bisection method for finding roots
% Author: Caleb Jacobs
% Date last modified: 09-09-2021
%
% a: left searching interval bound
% b: right searching interval bound
% f: function to compute the root of
% tol: stopping tolerance
% maxIts: maximum iterations allowed
% root: true root for error calculations
% return: error code
function ier = bisection(a, b, f, tol, maxIts, root)
    fprintf("Bisection Method\n")

    % Check if root is guaranteed in the interval
    if f(a) * f(b) > 0
        ier = 2;                        % Return 2 if solution not guaranteed
        return
    end
    
    % Begin iterations
    for i = 1:maxIts
        c = (b + a) / 2;                % Midpoint of interval
        
        fprintf('%3d: %.16e\n', i, abs(c - root) / abs(root))
        
        if f(c) == 0 || (b - c) < tol   % Check convergence criteria
            ier = 0;                    % Return 0 for success
            fprintf('Root = %.16e\n\n', c);
            return
        end
        
        if f(c) * f(a) >= 0             % Update searching interval
            a = c;
        else
            b = c;
        end
    end
    
    ier = 1;                            % Return 1 for no convergence
    return
end