%% Normalized Power Iteration %%
clear

% Preliminary settings
mits = 200; % max iterations
eps = 10^(-10);

% Define A
A = [6,2,1;
	2,3,1;
	1,1,1];

% Define x0
x(1:3,1) = [0;0;1];
data = [0, x(:,1)', 0];

% Begin iterations with a 
for k = 1:mits
	y(1:3,k) = A*x(:,k);
	my = max(abs(y(:,k)))
	x(1:3,k+1) = y(:,k)./my;
	
	% Formatting data
	data = [data; k, x(:,k+1)', my];
	
	% Convergence test
	if max(abs(x(:,k+1)-x(:,k))) < eps
		break
	end
end