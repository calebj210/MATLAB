% Homework 1 problem 2
% Author: Caleb Jacobs
% Date last modified 26-08-2021


% Polynomial defintions
p1 = @(x) (x - 2) .^ 9;
p2 = @(x) x.^9 - 18*x.^8 + 144*x.^7 - 672*x.^6 + 2016*x.^5 - 4032*x.^4 + 5376*x.^3 - 4608*x.^2 + 2304*x - 512;

x = [1.920 : 0.001 : 2.080];

f1 = figure;
plot(x, p1(x), 'LineWidth', 1.5)
axis square
title('p(x) = (x - 2)^9')
xlabel('x')
ylabel('p(x)')
saveas(f1, '2.ii.png')


f2 = figure;
plot(x, p2(x), 'LineWidth', 1.5)
axis square
title('p(x) = Expanded (x - 2)^9')
xlabel('x')
ylabel('p(x)')
saveas(f2, '2.i.png')