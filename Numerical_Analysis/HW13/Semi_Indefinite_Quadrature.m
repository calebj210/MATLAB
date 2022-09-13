%%
% Quadrature for computing the gamma function
% 
% Date last modified: 09-Dec-2021
% Author: Caleb Jacobs
format longe
gam(4, 1000) - gamma(4)


function val = gam(x, n)
    b = 20 * x;                     % Compute upper bound of integration 
    
    t = linspace(0, b, n);          % Create evaluation nodes
    h = t(2) - t(1);                % Find stepsize
    f = t .^ (x - 1) .* exp(-t);    % Compute function values at each node
    
    val = h * ((f(1) + f(n)) / 2 + sum(f(2 : n - 1)));
end