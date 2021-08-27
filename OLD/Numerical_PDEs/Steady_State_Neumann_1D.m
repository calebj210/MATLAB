%%% 1D Steady State Heat Equation with Neumann Solver %%%
clear
close all
format long
% f1 = figure;
% f2 = figure;

%% Parameters
% Number of points
pts = [20,40,80];

% Boundary conditions
a = 0;
alpha = 1;
b = 1;
beta = -1;
f = @(x) x.^2;

% Exact solution
trU = @(x) (x.^4+12*x-25)/12;

%% Solver for AU=F
errs = zeros(size(pts,2),1);
j = 1;
for k = pts
	M = k;
	h = 1/(M+1);
	
	% Construct A
	A = zeros(M+1,M+1);
	A(2:end,2:end) = diag(ones(M-1,1),1) + -2*eye(M) + diag(ones(M-1,1),-1);
	A(1,1:2) = [-h,h];
	A(2,1) = 1;
	A = A/(h^2);
	
	% Construct F
	F = zeros(M+1,1);
	F(1) = alpha + h*f(a)/2;
	for i = 2:M+1
		F(i) = f(a + (i-1)*h);
	end
	F(end) = F(end)-beta/(h^2);
	
	% Solve for U
	U = A\F;
	
	%% Analysis
	% Plot results
% 	figure(f1)
	t = a+h:h:b-h;
% 	subplot(3,1,j);
% 	plot(t,U(2:end),'k',t,trU(t),'o')
% 	title(strcat("Solution Using M = ", num2str(k)))
% 	legend('Approximate U','True U')
% 	xlabel('x')
% 	ylabel('u(x)')
% 	
% 	% Plot error vectors
% 	figure(f2)
% 	subplot(3,1,j)
% 	plot(t,U(2:end)-trU(t)')
% 	title(strcat('Error Using M = ', num2str(k)))
% 	xlabel('x')
% 	ylabel('U - True U')
	
	% Compute and store 2-norm of error
	errs(j) = norm(U(2:end)-trU(t)',inf);
	j = j+1;
end

%% Error table
clc
tbl = table(pts',errs,'VariableNames',["M","Norm_Error"])