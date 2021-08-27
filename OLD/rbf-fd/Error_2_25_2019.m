%%%%%%%%%%%%%%%%
%% My first try at RBF-FD %%
%%%%%%%%%%%%%%%%

%% Solving the heat equation over a line %%
clear;
close all;
relErr=zeros(5,91);

for k=3:7
for j=10:100
%% Setup %%
dt = .001; % Time step
endtime = 1; % End time
its = round(endtime/dt); % Number of iterations
N = j; % Approximate amount of nodes
N_BDR = 2; % Number of boundary nodes
n = k; % Number of neighbors
m=3;

%% RBF's Setup%%
%{
phi = @(x1,x2) exp(-m.*(x1-x2).^2);
phi_x = @(x1,x2) m.*(x1-x2)*((x1-x2).^2).^(m./2-1);
phi_xx = @(x1,x2) 2.*exp(-m.*(x1-x2).^2).*m.*(-1+2.*m.*(x1-x2).^2);
%}
phi = @(x1,x2) ((x1+(-1).*x2).^2).^((1/2).*m);
phi_x = @(x1,x2) 2.*exp(1).^((-1).*m.*(x1+(-1).*x2).^2).*m.*((-1).*x1+x2);
phi_xx = @(x1,x2) ((-1)+m).*m.*((x1+(-1).*x2).^2).^((1/2).*((-2)+m));


%% Constructing our Nodes %%
x = linspace(0,2,N)';
[idx,dist] = knnsearch(x,x,'k',n);

%% Boundaries %%
f = x.^2;
truef = 2.*ones(N,1);

%% Populating our Differentiation Matrix %%
Dlap = sparse(N);
for i = 1:N
	
	xn = x(idx(i,:),1);
	
	[X1 X2] = meshgrid(xn);
	
	A = phi(X1,X2);
	b = phi_xx(X1(1,:),X2(1,:));
	
	Dlap_local = b/A;
	
	Dlap(i,idx(i,:)) =  Dlap_local;
end

relErr(k-2,j-9)=norm(Dlap*f-truef)./norm(truef);
end
end
for k=1:5
	loglog(10:100,relErr(k,:))
	hold on
end
hold off








