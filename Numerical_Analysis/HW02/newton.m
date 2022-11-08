%%
% Newton's root finding method
% Author: Caleb Jacobs
% Date last modified: 07-09-2021
%
% x0: initial root guess
% f:  function to find root of
% df: derivative of function 
% maxIts: maximum iterations for convergence
% tol: convergence tolerance
% root: true root for error computation
% return: error code
function ier = newton(x0, f, df, maxIts, tol, root)
    fprintf("Newton's Method\n")
    
    % Begin running iterations
    for i = 1:maxIts
        denom = df(x0);             % Compute derivative of f
        
        if denom == 0               % Check if derivative is flat
            ier = 2;
            return
        end
        
        x1 = x0 - f(x0) / denom;    % Compute next iteration
        
        fprintf('%3d: %.16e\n', i, abs(x1 - root) / abs(root))
        
        if abs(x1 - x0) <= tol      % Check stopping criteria
            ier = 0;                % Return 0 for success
            fprintf('Root = %.16e\n\n', x1);
            return
        end
        
        x0 = x1;                    % Save current iteration
    end
    
    ier = 1;                        % Return 1 if method did not converge
    return
end