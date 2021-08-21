%%% Euler Forward and Backward for 1D %%%
clear
close all
format shorte

%% Settings
% Step size
step = [0.1,0.05,0.005];

% Method (Forward or Backward)
method = 'Backward';

% ODE
f = @(u,t) -100*u;
u0 = 4;
ti = 0;
tf = 2;

% True solution
trU = @(t) 4*exp(-100*t);

%% Begin iteration
j = 1;
errs = zeros(size(step))';
for k = step
its = (tf-ti)/k;
u = zeros(its+1,1);
u(1) = u0;
t = ti:k:tf;
for i = 1:its
	switch method
		case 'Forward'
			u(i+1) = u(i)+k*f(u(i),t(i));
		case 'Backward'
			u(i+1) = u(i)/(1+100*k);
	end
end

subplot(size(step,2),1,j)
plot(t,u,'O',t,trU(t)')
title(sprintf('%s Euler Method with k = %0.3f',method, k))
legend('Approximate Solution','True Solution')
xlabel('x')
ylabel('u(x)')

errs(j) = norm(u-trU(t)',1)/its;
j = j+1;
end

errs(2)/errs(3)