%%% Steady state solver using Chebychev polynomial interpolation
clear
close all
format long
f1 = figure;
f2 = figure;


%% Settings and boundary conditions
% Symbolic variables
syms xi;

% Number of grid points
pts = [10,20];

% Boundary values
a = -1;
b = 1;
alpha = 0;
beta = 1;
f = @(x) exp(x);

% True solution
trU = @(x) (-1/2).*exp(1).^(-1).*(1+(-1).*exp(1)+exp(1).^2+(-2).*exp(1).^(1+x)+(-1).*x+(-1).*exp(1).*x+exp(1).^2.*x);

% Compute solution using special Chebyshev grid
errsC = zeros(size(pts,2),1);
j = 1;
for k = pts
	M = k;
	
	row(xi) = chebyshevT(0:M+1,xi);
	A = zeros(M+2);
	A(1,:) = row(a);
	A(end,:) = row(b);
	
	x = -cos((1:M)*pi/(M+1));
	row2 = diff(row,2);
	for i = 1:M
		A(i+1,:) = row2(x(i));
	end
	
	F = zeros(M+2,1);
	F(1) = alpha;
	F(end) = beta;
	F(2:end-1) = f(x);
	
	C = A\F;
	
	calc = zeros(M,1);
	for i = 1:M
		calc(i) = dot(C,row(x(i)));
	end
	
	figure(f1)
	subplot(size(pts,2),2,2*j-1)
	plot(x,calc-trU(x)')
	title(strcat("Chebyshev Spacing Error: M = ", num2str(k)))
	xlabel('x')
	ylabel('Error')
	
	figure(f2)
	subplot(size(pts,2),2,2*j-1)
	plot(x,calc,'o',x,trU(x))
	title(strcat("Chebyshev Spacing Solution: M = ", num2str(k)))
	xlabel('x')
	ylabel('u(x)')
	
	errsC(j) = norm(calc-trU(x)')/sqrt(M);
	j = j+1;
end

% Compute solution using uniform grid
errsG = zeros(size(pts,2),1);
j = 1;
for k = pts
	M = k;
	h = (b-a)/(M+1);
	
	row(xi) = chebyshevT(0:M+1,xi);
	A = zeros(M+2);
	A(1,:) = row(a);
	A(end,:) = row(b);
	
	x = a+h:h:b-h;
	row2 = diff(row,2);
	for i = 1:M
		A(i+1,:) = row2(x(i));
	end
	
	F = zeros(M+2,1);
	F(1) = alpha;
	F(end) = beta;
	F(2:end-1) = f(x);
	
	C = A\F;
	
	calc = zeros(M,1);
	for i = 1:M
		calc(i) = dot(C,row(x(i)));
	end
	
	figure(f1)
	subplot(size(pts,2),2,2*j)
	plot(x,calc-trU(x)')
	title(strcat("Uniform Spacing Error: M = ", num2str(k)))
	xlabel('x')
	ylabel('Error')
	
	figure(f2)
	subplot(size(pts,2),2,2*j)
	plot(x,calc,'o',x,trU(x))
	title(strcat("Uniform Spacing Solution: M = ", num2str(k)))
	xlabel('x')
	ylabel('u(x)')
	
	
	errsG(j) = norm(calc-trU(x)')/sqrt(M);
	j = j+1;
end

% Order 2 FD solution
errsFD = zeros(size(pts,2),1);
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
	t = a+h:h:b-h;
	
	% Compute and store 2-norm of error
	errsFD(j) = norm(U-trU(t)')/sqrt(M);
	j = j+1;
end

%% Error table
clc
tbl = table(pts',errsC,errsG,errsFD,'VariableNames',["M","Chebyshev_Spacing","Uniform_Grid","Second_Ord_FD"])