%%% Burgers equation solver using finite volume method with Lax-Friedrichs
%%% flux

clear
close all
format long
clc

stls = strings(1,5);
stls(1) = '.k';stls(2) = ':k';stls(3) = '-.k';stls(4) = '--k';stls(5) = '-k';

grd = [20,40,80,160,320,640,1280,2560];
grd = 20*2.^[0:4];
errs = zeros(size(grd,2),1);
for ii = 1:size(grd,2)

%% Settings
% Cell average setup
N = grd(ii);
a = 0;
b = 2*pi;

% Final time
T = 1;

% CFL condition
omega = 0.5;

% Initial conditions 
u0 = @(x) 1+ 0.5*sin(x);
u0p = @(x) 0.5*cos(x);

% Average value integral
uVol = @(x1,x2) 1 + (cos(x1) - cos(x2))/(2*(x2-x1));

%% True solution settings
% Max iterations to use
maxIt = 1000;

% Tolerance
epsilon = 10^(-13);

% Initial conditons and derivatives
f = @(x,x0,t) x0 + u0(x0)*t - x;
ff = @(x0,t) 1 + u0p(x0)*t;

%% Setup
h = (b-a)/N;
x = linspace(a+h/2,b-h/2,N);

% Compute U0 averages
U = zeros(N,1);
for i = 1:N
	U(i) = uVol(x(i)-h/2,x(i)+h/2);
end
% U = u0(x(1:N));

% Initialize time variable
t = 0;

%% Begin solving
% plot(x,U)
% pause(0.1)
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
	
% 	plot(x,U)
% 	ylim([0.5,1.5])
% 	pause(0.001)
end

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

errs(ii) = norm(U-trU,1)/N;

% plot(x,U,stls(ii),'DisplayName',strcat('N = ',num2str(grd(ii))))
% legend('Location','northwest')
% title('Solution at T = 1 using various N')
% xlabel('x')
% ylabel('u(x)')
% hold on
end

% plot(x,trU,'-d','DisplayName','True Solution')

cnvgnc = NaN(size(grd,2),1);
for i = 2:size(grd,2)
	cnvgnc(i) = log2(errs(i-1)/errs(i));
end

table(grd',errs,cnvgnc,'VariableNames',["N","Norm(Error)","Convergence Rate"])