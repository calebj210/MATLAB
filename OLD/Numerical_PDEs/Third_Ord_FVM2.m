%%% Third order FVM Solver %%%

clear
close all
clc
format long

% Max time
T = 3;

% Boundaries
a = 0;
b = 2*pi;

%% Number of grid points to use (can add multiple points to get convergence rates and error comparisons)
Grids = [160];

% Preallocate space for errors
errs = Grids;

% Set different grid sizes
for s = 1:size(Grids,2)
% Spatial grid points
N = Grids(s);

% Compute delta x
h = (b-a)/N;

% Initialization of grid
x = NaN(N+1,1);
for i = 1:N+1
	x(i) = a + (i-1)*h;
end

%% Step function initial condition
u = NaN(N+4,1);
for i = 3:N+2
	j = i - 2;
	if x(j) < pi
		u(i) = 0;
	else
		u(i) = 1;
	end
end

% Setup periodic boundary condtions
u(1) = u(N+1);
u(2) = u(N+2);
u(end-1) = u(3);
u(end) = u(4);

% CFL condition setup
cfl = 0.6;
k = cfl*h;
M = floor(T/k)+1;
k = T/M;

%% TVD RK3
% Preaclloacte space for u_i+1/2 and the temp u
uh = u;
utmp = u;

% Initial error plot at each node
% plot(x(1:end-1),u(3:N+2) - -(cos(x(2:end)) - cos(x(1:end-1)))/h)
% pause(0.5)

% Begin time stepping
for t = 1:M
	for i = 3:N+2
		
		% Adaptive stencil choice
		l = 0;
		if abs(u(i) - u(i-1)) <= abs(u(i+1) - u(i))
			l = l-1;
		end
		if abs(u(i+l-1) - 2*u(i+l) + u(i+l+1)) <= abs(u(i+l) - 2*u(i+l+1) + u(i+l+2))
			l = l-1;
		end
		
		% Compute polynomial reconstruction
		switch l
			case -2
				%u^(1)_(i+1/2)
				uh(i) = (1/3)*u(i-2) - (7/6)*u(i-1) + (11/6)*u(i);
			case -1
				%u^(2)_(i+1/2)
				uh(i) = -(1/6)*u(i-1) + (5/6)*u(i) + (1/3)*u(i+1); 
			case 0
				%u^(3)_(i+1/2)
				uh(i) = (1/3)*u(i) + (5/6)*u(i+1) - (1/6)*u(i+2);
		end

	end
	uh(2) = uh(N+2);


	%%%% Implicit 3rd order Runge-Kutta step 1 %%%%
	lambda = k/h;
	for i = 3:N+2
		utmp(i) = u(i) - lambda*(uh(i) - uh(i-1));
	end
	
	
	% Intermediate boundary condtions
	utmp(1) = utmp(N+1);
	utmp(2) = utmp(N+2);
	utmp(end-1) = utmp(3);
	utmp(end) = utmp(4);
	
	% Second reconstruction
	for i = 3:N+2
		
		% Adaptive stencil choice
		l = 0;
		if abs(utmp(i) - utmp(i-1)) <= abs(utmp(i+1) - utmp(i))
			l = l-1;
		end
		if abs(utmp(i+l-1) - 2*utmp(i+l) + utmp(i+l+1)) <= abs(utmp(i+l) - 2*utmp(i+l+1) + utmp(i+l+2))
			l = l-1;
		end
		
		% Compute polynomial reconstruction
		switch l
			case -2
				%u^(1)_(i+1/2)
				uh(i) = (1/3)*utmp(i-2) - (7/6)*utmp(i-1) + (11/6)*utmp(i);
			case -1
				%u^(2)_(i+1/2)
				uh(i) = -(1/6)*utmp(i-1) + (5/6)*utmp(i) + (1/3)*utmp(i+1); 
			case 0
				%u^(3)_(i+1/2)
				uh(i) = (1/3)*utmp(i) + (5/6)*utmp(i+1) - (1/6)*utmp(i+2);
		end

	end
	uh(2) = uh(N+2);
	
	%%%% Implicit 3rd order Runge-Kutta step 2 %%%%
	for i = 3:N+2
		utmp(i) = (3/4)*u(i) + (1/4)*(utmp(i) - lambda*(uh(i) - uh(i-1)));
	end
	
	% Intermediate boundary condtions
	utmp(1) = utmp(N+1);
	utmp(2) = utmp(N+2);
	utmp(end-1) = utmp(3);
	utmp(end) = utmp(4);
	
	% Third reconstruction
	for i = 3:N+2
		
		% Adaptive stencil choice
		l = 0;
		if abs(utmp(i) - utmp(i-1)) <= abs(utmp(i+1) - utmp(i))
			l = l-1;
		end
		if abs(utmp(i+l-1) - 2*utmp(i+l) + utmp(i+l+1)) <= abs(utmp(i+l) - 2*utmp(i+l+1) + utmp(i+l+2))
			l = l-1;
		end
		
		% Compute polynomial reconstruction
		switch l
			case -2
				%u^(1)_(i+1/2)
				uh(i) = (1/3)*utmp(i-2) - (7/6)*utmp(i-1) + (11/6)*utmp(i);
			case -1
				%u^(2)_(i+1/2)
				uh(i) = -(1/6)*utmp(i-1) + (5/6)*utmp(i) + (1/3)*utmp(i+1); 
			case 0
				%u^(3)_(i+1/2)
				uh(i) = (1/3)*utmp(i) + (5/6)*utmp(i+1) - (1/6)*utmp(i+2);
		end

	end
	uh(2) = uh(N+2);
	
	%%%% Implicit 3rd order Runge-Kutta step 3 %%%%
	for i = 3:N+2
		u(i) = (1/3)*u(i) + (2/3)*(utmp(i) - lambda*(uh(i) - uh(i-1)));
	end
	
	% Final boundary conditions
	u(1) = u(N+1);
	u(2) = u(N+2);
	u(end-1) = u(3);
	u(end) = u(4);
	
	
%% Error at each node over time (Just uncomment)
% 	pause(0.01)
% 	plot(x(1:end-1),u(3:N+2) - -(cos(x(2:end) - t*k) - cos(x(1:end-1) - t*k))/h)
end

%% True solution calculation
trU = heaviside(x(1:N) - pi - T) + heaviside(-(x(1:N) - T));

%% Plot of solution and exact solution
plot(x(1:N),u(3:N+2),'o', x(1:N),trU)
title('ENO Scheme')
xlabel('x')
ylabel('u(x,T)')
legend('Numerical Solution','Exact Solution')

%% Error calculation
errs(s) = norm(u(3:N+2)-trU,1)*h;
end

%% Convergence rate calculation
cnvgs = NaN(size(Grids,2),1);
for i = 2:size(Grids,2)
	cnvgs(i) = log2(errs(i-1)/errs(i));
end

%% Display Table
table(Grids',errs',cnvgs,'VariableNames',["N","Norm(Error)","Convergence Rate"])