%%% 1D Burger's Equation Solver (Upwinding)%%%
clear
close all
format long

%% Setting
% Time step and final time
k = 0.001;
tf = 1;

% Grid points
M = 1000;

% Initial condition
f0 = @(x) 1/2 + sin(x);

%% Begin solving
T = tf/k;
h = 2*pi/(M+1);
u0 = f0(0:h:2*pi);
u = zeros(size(u0));
U = zeros(T,M+2);
U(1,:) = u0;

% plot(0:h:2*pi,u0)
% title(0)
% pause(2)
for t = 1:T+1
	
	% Solve first node
	if u0(1) >= 0
		u(1) = u0(1) - (k/(2*h))*(u0(1)^2-u0(M+2)^2);
	else
		u(1) = u0(1) - (k/(2*h))*(u0(2)^2-u0(1)^2);
	end
	
	% Solve intermediate nodes
	for j = 2:M+1
		if u0(j) >= 0
			u(j) = u0(j) - (k/(2*h))*(u0(j)^2-u0(j-1)^2);
		else
			u(j) = u0(j) - (k/(2*h))*(u0(j+1)^2-u0(j)^2);
		end
	end
	
	% Solve end node
	if u0(M+2) >= 0
		u(M+2) = u0(M+2) - (k/(2*h))*(u0(M+2)^2-u0(M+1)^2);
	else
		u(M+2) = u0(M+2) - (k/(2*h))*(u0(1)^2-u0(M+2)^2);
	end
	
	plot(0:h:2*pi,u)
	title((t)*k)
	pause(0.001)
	u0 = u;
	U(t+1,:) = u;
end

% [X,Y] = meshgrid(0:h:2*pi,0:k:tf+k);
% surf(X,Y,U,'EdgeColor','none')
% xlim([0 2*pi])
% ylim([0 tf])
% zlim([-0.5 1.5])
% xlabel('x')
% ylabel('t')
% zlabel('u(x,t)')
% colorbar
% colormap gray;
% title('Approximate Solution')