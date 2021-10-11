%%
% 2D polynomial interpolation
% Author: Caleb Jacobs
% Date last modified: 11-10-2021

function Problem_2
    %% Settings
    format longE

    %% Setup
    f  = @(x) exp(x(:,1)) .* sin(x(:,2));       % Function to interpolate
    xi = [0, 0; ...
          0, 2; ...
          1, 0; ...
          1, 2; ...
          2, 1; ...
          2, 3];                                % Nodes
    n = length(xi);                             % Number of nodes
    
    %% Create system of equations
    % Create LHS matrix
    A = [ones(n,1), xi(:,1), xi(:,2), xi(:,1) .* xi(:,2), xi(:,1).^2, xi(:,2).^2];
    b = f(xi);                                  % RHS function vector
    
    %% Find polynomial coefficients
    c = A \ b
    
    %% Setup for plotting interpolation
    [X, Y] = meshgrid(linspace(-1, 3, 25), linspace(-1, 3, 25));
    F = @(x,y) exp(x) .* sin(y);                % Original function
    
    P = @(x,y) c(1) + c(2)*x + c(3)*y + c(4)*x.*y + c(5)*x.^2 + c(6)*y.^2;
    
    %% Plot surfaces
    figure(1), clf
    surf(X, Y, P(X,Y))
    title('Polynomial Interpolant of f')
    zlim([-20.1,20.1])
    
    figure(2), clf
    surf(X, Y, F(X,Y))
    title('True f')
    zlim([-20.1,20.1])
    
end