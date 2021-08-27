%%% Newton's Method Equation Solver %%%
clear
close all
format long

%% Settings
% Tolerance
epsilon = 10^(-7);
% Initial guesses
% init = [-2,-0.2,0.2,2];
init = [-2,-0.2,0.2,2];

% Max iterations allowed
maxIt = 100;

% Function definition
f = @(x) x^3 - sin(x);
ff = @(x) 3*x^2 - cos(x);


j = 1;
calc = zeros(size(init,2),1);
itsUsed = calc;
for k = init
	x0 = k;
	
	for i = 0:maxIt
		if abs(f(x0)) <= epsilon
			break
		end
		
		x0 = x0 - f(x0)/ff(x0)
	end
	
	calc(j) = x0;
	itsUsed(j) = i;
	j = j+1;
end

clc
table(init',calc,itsUsed,'VariableNames',["x0","Solution_(x)","Iterations_Used"])