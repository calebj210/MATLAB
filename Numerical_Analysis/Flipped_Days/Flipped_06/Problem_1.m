%% 
% Compare orthogonal legendre polynomial bases
% 
% Author: Caleb Jacobs
% Date last modified: 07-11-2021

%% Parameters
n = 100000;

%% Compute polynomials
x = linspace(-1, 1, 3000)';

tic
p = LP(x, n);
toc

%% Plot all polynomials
figure(1)
for i = n-5 : n + 1
    plot(x, p(:, i), 'LineWidth', 2)
    hold on
end
xlim([-1, 1])
ylim([-2, 2])
hold off

%% Plot last polynomial
figure(2)
plot(x, p(:, end), 'LineWidth', 2)
xlim([-1, 1])
ylim([-2, 2])

%% Evaluate Legendre Polynomials
function y = LP(x, N)
    y = zeros(length(x), N + 1);        % Initialize polynomial
    y(:, 1) = 1;                        % Compute 0th polynomial
    y(:, 2) = 3 * x / 2;                % Compute 1st polynomial
    
    % Recursively compute the other polynomials
    for i = 3 : N + 1
        n = i - 1;
        y(:, i) = ((2*n + 1) / (n + 1)) * x .* y(:, i - 1) - (n / (n + 1)) * y(:, i - 2);
    end
end