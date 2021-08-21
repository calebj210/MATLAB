%%% Conjugate Gradient Iteration for 1D %%%
clear
close all
format long

%% Settings
% Tolerance
epsilon = 10^(-8);

% Max iterations
maxIter = 500;

% Number of nodes
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

% Construct F
F = zeros(m,1);
F(1) = f(a+h) - alpha/(h^2);
F(m) = f(b-h) - beta/(h^2);
for i = 2:m-1
	F(i) = f(a + i*h);
end

% Begin iterations
u0 = zeros(size(F));
r0 = F - A*u0;
p0 = r0;
for i = 0:maxIter
	if abs(A*u0-F) <= epsilon
		break
	end
	
	alpa = (r0'*p0)/(p0'*A*p0);
	u0 = u0 + alpa*p0;
	r0 = r0 - alpa*A*p0;
	bta = -(r0'*A*p0)/(p0'*A*p0);
	p0 = r0 + bta*p0;
end

x = linspace(a+h,b-h,size(u0,1));
subplot(size(nodes,2),1,k)
plot(x,u0,'O',x,trU(x))
title(sprintf('%s Iteration with M = %i','Conjugate-Gradient',m))
xlabel('x')
ylabel('u(x)')

its(k) = i;
k = k+1;
end

its