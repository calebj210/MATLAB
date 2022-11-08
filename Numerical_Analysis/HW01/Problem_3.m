% Homework 1 problem 3
% Author: Caleb Jacobs
% Date last modified 26-08-2021

% Function defintions
f1 = @(x, delta) cos(x + delta) - cos(x);
f2 = @(x, delta) -2 * sin((2*x + delta) / 2) .* sin(delta / 2);
f3 = @(x, delta) -delta * sin(x) - (delta .^ 2) / 2 * cos(x);

% Variable defintions
X = logspace(-16,10,10000);
% X = linspace(0, 10000, 1000);
startColor = [0.5, 0.5, 0.5];
D = logspace(-16, 0, 17);

close all

% Plot differences
figure
for i = 1:length(D)
    loglog(X, abs(f1(X, D(i)) - f2(X, D(i))), 'LineWidth', 1.5, 'Color', startColor + (i - 1) / 32)
    hold on
end
title('Original vs My Expression')
ylim([1e-40 1e0])
xlabel('x')
ylabel('Abs error')
hold off

figure
for i = 1:length(D)
    loglog(X, abs(f3(X, D(i)) - f2(X, D(i))), 'LineWidth', 1.5, 'Color', startColor + (i - 1) / 32)
    hold on
end
title('My Expression vs Taylor')
ylim([1e-40 1e0])
xlabel('x')
ylabel('Abs error')
hold off

figure
for i = 1:length(D)
    loglog(X, abs(f3(X, D(i)) - f1(X, D(i))), 'LineWidth', 1.5, 'Color', startColor + (i - 1) / 32)
    hold on
end
title('Original vs Taylor')
ylim([1e-40 1e0])
xlabel('x')
ylabel('Abs error')
hold off