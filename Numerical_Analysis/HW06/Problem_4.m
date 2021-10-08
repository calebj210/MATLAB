%%
% Homework 6, problem 4 code
% Newton like iteration
%
% Author: Caleb Jacobs
% Date last modified: 07-10-2021

clear
close all
format long

%% Setup
f1 = @(x) 3*x(1).^2 + 4*x(2).^2 - 1;
f2 = @(x) x(2).^3 - 8*x(1).^3 - 1;
g  = @(x) x - [0.016, -0.17; 0.52, -0.26] * [f1(x); f2(x)];

%% Begin fixed point iterations
x0 = [-0.5; 0.25];
for i = 1:100
    x0 = g(x0);
    if abs(f1(x0) - f2(x0)) < 1e-7
        break
    end
end

%% Display information
i
x0
f1(x0)
f2(x0)