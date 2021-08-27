%%% 1D Steady State Heat Equation with Dirchlet BC Solver %%%
clear
close all
format long
f1 = figure;
% f2 = figure;

%% Parameters
% Number of points
pts = [10,20,40];

% Boundary conditions
% a = 0;
% alpha = 0;
% b = 1;
% beta = 0;
% f = @(x) x;
% 
% % Exact solution
% trU = @(x) (x.^3-x)/6;

a = -1;
b = 1;
alpha = 0;
beta = 1;
f = @(x) exp(x);

% True solution
trU = @(x) (-1/2).*exp(1).^(-1).*(1+(-1).*exp(1)+exp(1).^2+(-2).*exp(1).^(1+x)+(-1).*x+(-1).*exp(1).*x+exp(1).^2.*x);


%% Solver for AU=F
errs = zeros(size(pts,2),1);
j = 1;
for k = pts
	M = k;
	h = (b-a)/(M+1);
	
	% Construct A
	A = diag(ones(M-1,1),1) + -2*eye(M) + diag(ones(M-1,1),-1);
	A = A/(h^2);
	
	% Construct F
	F = zeros(M,1);
	F(1) = f(a+h) - alpha/(h^2);
	F(M) = f(b-h) - beta/(h^2);
	for i = 2:M-1
		F(i) = f(a + i*h);
	end
	
	% Solve for U
	U = A\F;
	
	% Plot results
	figure(f1)
	t = a+h:h:b-h;
	subplot(3,1,j);
	plot(t,U,'k',t,trU(t),'o')
	title(strcat("Finite Difference Solution Using M = ", num2str(k)))
	legend('Approximate U','True U')
	xlabel('x')
	ylabel('u(x)')
% 	
% 	% Plot error vectors
% 	figure(f2)
% 	subplot(3,1,j)
% 	plot(t,U-trU(t)')
% 	title(strcat('Error Using M = ', num2str(k)))
% 	xlabel('x')
% 	ylabel('U - True U')
	
	% Compute and store 2-norm of error
	errs(j) = norm(U-trU(t)',1)/M;
	j = j+1;
end

clc
tbl = table(pts',errs,[NaN;log2(errs(1)/errs(2));log2(errs(2)/errs(3))],'VariableNames',["M","Norm_Error","Convgence_Rate"])