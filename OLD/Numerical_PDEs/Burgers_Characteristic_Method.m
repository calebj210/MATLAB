%%% Numerical Burgers Equation Solver Using Characteristics
clear
close all
format long

%% Setting
% Tolerance and max iterations
epsilon = 10^(-7);
maxIt = 10000;

% Time step and final time
k = 0.001;
tf = 1;

% Grid points
a = 0;
b = 2;
M = 40;

% Initial condition
u0 = @(x) 1 + tanh(x);
u0p = @(x) sech(x).^2;

% Characteristic equation
f = @(x,x0,t) x0 + u0(x0)*t - x;
ff = @(x0,t) 1 + u0p(x0)*t;

%% Solver
X = linspace(a,b,M);
t = linspace(0,tf,tf/k);
U = zeros(1,M);

% Full solution
totU = zeros(tf/k,M);

% Begin solving
% Time stepper
for k = 1:(tf/k)
	% Solve at each x
	for i = 1:M
		x = X(i);
		
		% Initial guess for newtons method
		x0 = x;

		% Solve characteristic equation using newtons method.
		for j = 0:maxIt
			if abs(f(x,x0,t(k))) <= epsilon
				break
			end
			
			x0 = x0 - f(x,x0,t(k))/ff(x0,t(k));
		end
		
		U(i) = u0(x0);
	end 
	
	plot(X,U)
	pause(0.001)
	totU(k,:) = U;
end
plot(X,U)

[X,T] = meshgrid(X,t);
subplot(2,1,1)
surf(X,T,totU,'EdgeColor','none')
xlim([0 b])
ylim([0 tf])
zlim([-0.5 1.5])
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
c = colorbar;
c.Label.String = 'u(x,t)';
title('Approximate Solution Surface Plot')
subplot(2,1,2)
contourf(X,T,totU,20)
xlim([0 b])
ylim([0 tf])
xlabel('x')
ylabel('t')
c = colorbar
c.Label.String = 'u(x,t)';
colormap gray;
title('Approximate Solution Contour Map')