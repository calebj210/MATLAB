%%% Gauss-Seidel and Jacobi Iteration Methods in 1D %%%
%%% I also implemented an SOR iteration scheme %%%
clear
close all
format long


%% Settings
% Tolerance
epsilon = 10^(-8);

% Method to use (Jacobi or Gauss-Seidel or SOR)
method = 'Jacobi';

% Max iterations
maxIter = 100000;

nodes = [20,40,80];

% Problem
a = -1;
alpha = 0;
b = 1;
beta = 1;

f = @(x) exp(x);

% True solution
trU = @(x) (-1/2).*exp(1).^(-1).*(1+(-1).*exp(1)+exp(1).^2+(-2).*exp(1).^(1+x)+(-1).*x+(-1).*exp(1).*x+exp(1).^2.*x);

%% Begin solving
k = 1;
its = zeros(size(nodes))';
for m = nodes
% Construct A
h = (b-a)/(m+1);
A = diag(ones(m-1,1),1) + -2*eye(m) + diag(ones(m-1,1),-1);
A = A/(h^2);

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
switch method
	case 'Jacobi'
		M = D;
		N = -L - U;
	case 'Gauss-Seidel'
		M = D + L;
		N = -U;
	case 'SOR'
		% Relaxation Coefficient
        omega = 1.75;
		M = 1/omega * (D + omega*L);
		N = 1/omega * ((1-omega)*D - omega*U);
end

% Construct F
F = zeros(m,1);
F(1) = f(a+h) - alpha/(h^2);
F(m) = f(b-h) - beta/(h^2);
for i = 2:m-1
	F(i) = f(a + i*h);
end

% Iterate
u0 = zeros(size(F));
for i = 0:maxIter
	if abs(A*u0-F) <= epsilon
		break
	end
	
	u0 = M\(N*u0 + F);
end

x = linspace(a+h,b-h,size(u0,1));
subplot(size(nodes,2),1,k)
plot(x,u0,'O',x,trU(x))
title(sprintf('%s Iteration with M = %i',method,m))
xlabel('x')
ylabel('u(x)')

its(k) = i;
k = k+1;
end

its