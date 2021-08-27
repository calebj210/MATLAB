%%% Burger's Equation Solver Using Finite Difference with Lax-Friedrich's
%%% flux with periodic boundary conditions

clear
close all
format long
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings %%
% Final time
T = 2;

% Spatial grid discretization
N = 160;

% CFL condition setting
omega = 0.9;

% Initial conditions and BC
a = -2;
b = 2;
l = -0.5;
r = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup %%
% Construct grid
h = (b-a)/N;
x = linspace(a+h/2,b-h/2,N);

% Initialize time variable
t = 0;

% Compute initial conditions
U = zeros(N,1);
for i = 1:N
	if x(i) < 0
		U(i) = l;
	else
		U(i) = r;
	end
end

% Begin time iteration
% plot(x,U)
% pause(1)
while t < T
	% Store previous step
	U0 = U;
	
	% Compute alpha
	alpha = max(U0);
	
	% Modify time step
	k = (h/alpha)*omega;
	k = T/(floor(T/k)+1);
	
	% Compute left endpoint
	if U0(1) >= 0
		uMinus = l;
	else
		uMinus = U0(1);
	end
	fplus = (1/2)*((U0(1)^2)/2 + (U0(2)^2)/2 - alpha*(U0(2) - U0(1)));
	fminus = (1/2)*((uMinus^2)/2 + (U0(1)^2)/2 - alpha*(U0(1) - uMinus));
	U(1) = U0(1) - (k/h)*(fplus - fminus);
	
	% Copmute intermediate points
	for j = 2:N-1
		fplus = (1/2)*((U0(j)^2)/2 + (U0(j+1)^2)/2 - alpha*(U0(j+1) - U0(j)));
		fminus = (1/2)*((U0(j-1)^2)/2 + (U0(j)^2)/2 - alpha*(U0(j) - U0(j-1)));
		U(j) = U0(j) - (k/h)*(fplus - fminus);
	end
	
	% Compute right endpoint
	if U0(1) <= 0
		uPlus = r;
	else
		uPlus = U0(N);
	end
	fplus = (1/2)*((U0(N)^2)/2 + (uPlus^2)/2 - alpha*(uPlus - U0(N)));
	fminus = (1/2)*((U0(N-1)^2)/2 + (U0(N)^2)/2 - alpha*(U0(N) - U0(N-1)));
	U(N) = U0(N) - (k/h)*(fplus - fminus);
	
	t = t + k;
	
	plot(x,U)
	ylim([-0.7,1.2])
	pause(0.01)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

