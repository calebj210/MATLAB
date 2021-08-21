clear
close all
format long
clc

tmp2 = zeros(3,1);
n = [160];
for s = 1:1
%%% Burger's Equation Solver Using Finite Difference with Lax-Friedrich's
%%% flux with periodic boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings %%
% Final time
T = 3;

% Spatial grid discretization
N = n(s);

% CFL condition setting
omega = 0.9;

% Initial conditions and BC
u0 = @(x) 1 + sin(x)/2;
a = 0;
b = 2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% True solution settings %%
% tolerance
epsilon = 10^(-13);

% Max iterations
maxIt = 10000;

% Initial condition derivative
u0p = @(x) cos(x);

% Characteristic equation
f = @(x,x0,t) x0 + u0(x0)*t - x;
ff = @(x0,t) 1 + u0p(x0)*t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup %%
% Construct grid
h = (b-a)/N;
x = linspace(a+h/2,b-h/2,N);

% Initialize time variable
t = 0;

% Compute initial conditions
U = u0(x);

% Begin time iteration
plot(x,U)
pause(1)
while t < T
	% Store previous step
	U0 = U;
	
	% Compute alpha
	alpha = max(U0);
	
	% Modify time step
	k = (h/alpha)*omega;
	k = T/(floor(T/k)+1);
	
	% Compute left endpoint
	fplus = (1/2)*((U0(1)^2)/2 + (U0(2)^2)/2 - alpha*(U0(2) - U0(1)));
	fminus = (1/2)*((U0(N)^2)/2 + (U0(1)^2)/2 - alpha*(U0(1) - U0(N)));
	U(1) = U0(1) - (k/h)*(fplus - fminus);
	
	% Copmute intermediate points
	for j = 2:N-1
		fplus = (1/2)*((U0(j)^2)/2 + (U0(j+1)^2)/2 - alpha*(U0(j+1) - U0(j)));
		fminus = (1/2)*((U0(j-1)^2)/2 + (U0(j)^2)/2 - alpha*(U0(j) - U0(j-1)));
		U(j) = U0(j) - (k/h)*(fplus - fminus);
	end
	
	% Compute right endpoint
	fplus = (1/2)*((U0(N)^2)/2 + (U0(1)^2)/2 - alpha*(U0(1) - U0(N)));
	fminus = (1/2)*((U0(N-1)^2)/2 + (U0(N)^2)/2 - alpha*(U0(N) - U0(N-1)));
	U(N) = U0(N) - (k/h)*(fplus - fminus);
	
	t = t + k;
	
	plot(x,U)
	ylim([0.5,1.5])
	pause(0.1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(x,U)
% legend('T = 0.5','T = 1','T = 1.5')
% hold on

%% Compute true solution %%
% True solution
trU = U;
for i = 1:N
		xi = x(i);
		
		% Initial guess for newtons method
		x0 = xi;

		% Solve characteristic equation using newtons method.
		for j = 0:maxIt
			if abs(f(xi,x0,T)) <= epsilon
				break
			end
			
			x0 = x0 - f(xi,x0,T)/ff(x0,T);
		end
		
		trU(i) = u0(x0);
end

% tmp2(s) = max(norm(U-trU,1))/N;
end

% table(n',tmp2,[NaN;log2(tmp2(2)/tmp2(1));log2(tmp2(3)/tmp2(2))],'VariableNames',["N","Norm(Error)","Convergence Rate"])