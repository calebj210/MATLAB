%%% Burgers equation solver using finite volume method with Lax-Friedrichs
%%% flux

clear
close all
format long
clc

stngs = strings(6,1);
stngs(1) = ':k';
stngs(2) = '-.k';
stngs(3) = '--k';
stngs(4) = '-k';

tms = [0,1,2,2.5];
for ii = 1:size(tms,2)
%% Settings
% Cell average setup
N = 160;
a = 0;
b = 2*pi;

% Final time
T = tms(ii);
k = 0.001;

% Initial conditions 
u0 = @(x) 1+ 0.5*sin(x);

% Average value integral
uVol = @(x1,x2) 1 + (cos(x1) - cos(x2))/(2*(x2-x1));

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
	
	% Compute left endpoint
	fplus = (((U0(1)+U0(2))/2)^2)/2;
	fminus = (((U0(N)+U0(1))/2)^2)/2;
	U(1) = U0(1) - (k/h)*(fplus - fminus);
	
	% Copmute intermediate points
	for j = 2:N-1
		fplus = (((U0(j)+U0(j+1))/2)^2)/2;
		fminus = (((U0(j)+U0(j-1))/2)^2)/2;
		U(j) = U0(j) - (k/h)*(fplus - fminus);
	end
	
	% Compute right endpoint
	fplus = (((U0(N)+U0(1))/2)^2)/2;
	fminus = (((U0(N)+U0(N-1))/2)^2)/2;
	U(N) = U0(N) - (k/h)*(fplus - fminus);
	
	t = t + k;
	
% 	plot(x,U)
% 	ylim([0.5,1.5])
% 	pause(0.001)
end

plot(x,U,stngs(ii),'DisplayName',['T = ',num2str(T)])
legend
xlabel('x')
ylabel('u(x)')
title('2nd Order(Method 1) Solution Over Time')
hold on
end