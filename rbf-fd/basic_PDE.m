%%%%%%%%%%%%%%%%
%% My first try at RBF-FD %%
%%%%%%%%%%%%%%%%

%% Solving the heat equation over a line %%
clear;
close all;

%% Setup %%
dt = .001; % Time step
endtime = 1; % End time
its = round(endtime/dt); % Number of iterations
N = 1000; % Approximate amount of nodes
N_BDR = 2; % Number of boundary nodes
n = 3; % Number of neighbors
m=1;

%% RBF's Setup%%
phi = @(x1,x2) exp(1).^((-1).*m.*(x1+(-1).*x2).^2);
phi_x = @(x1,x2) 2.*exp(1).^((-1).*m.*(x1+(-1).*x2).^2).*m.*((-1).*x1+x2);
phi_xx = @(x1,x2) 2.*exp(1).^((-1).*m.*(x1+(-1).*x2).^2).*m.*((-1)+2.*m.*(x1+(-1).*x2).^2);
%{
phi = @(x1,x2) exp(-m.*(x1-x2).^2);
phi_x = @(x1,x2) m.*(x1-x2)*((x1-x2).^2).^(m./2-1);
phi_xx = @(x1,x2) 2.*exp(-m.*(x1-x2).^2).*m.*(-1+2.*m.*(x1-x2).^2);
%}

%% Constructing our Nodes %%
x = linspace(0,2,N)';
[idx,dist] = knnsearch(x,x,'k',n);

%% Boundaries %%
f = x.^2;

%% Populating our Differentiation Matrix %%
Dlap = sparse(N);
for i = 1:N
	
	xn = x(idx(i,:),1);
	
	[X1,X2] = meshgrid(xn);
	
	A = phi(X1,X2);
	b = phi_xx(X1(1,:),X2(1,:));
	
	Dlap_local = b/A;
	
	Dlap(i,idx(i,:)) =  Dlap_local;
end

D = speye(N)-dt*Dlap;
D(1,:)=0;
D(N,:)=0;
D(1,1)=1;
D(N,N)=1;

%% Time Solver %%
u=f;
for i = 1:its
	u(1)=0;
	u(N)=0;
	
	u = D\u;
	
		figure(1)
	hold off
	plot(x,u)
	axis equal;
	grid on;
	axis([0,2,0,4])
	drawnow
end
