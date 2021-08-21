%% Raleigh Quotient Iteration %%
clear

% Preliminary settings
mits = 200; % max iterations
eps = 10^(-10);
B = 2; % Shift

% Define A
A = [6,2,1;
	2,3,1;
	1,1,1];

% Define x0
x(1:3,1) = [0;0;1];

% Calculate first Rayleigh Quotient
sig(1) = x(:,1)'*A*x(:,1)./(x(:,1)'*x(:,1));
data = [0, x(:,1)', 0, sig(1)];


% Begin iterations with a 
for k = 1:mits
	
	y(1:3,k) = (A-sig(k)*eye(3))\x(:,k);
	my = max(abs(y(:,k)));
	x(1:3,k+1) = y(:,k)./my;
	
	% Find next Rayleigh Quotient
	sig(k+1) = x(:,k+1)'*A*x(:,k+1)./(x(:,k+1)'*x(:,k+1));
	
	% Formatting results
	data = [data; k, x(:,k+1)', my, sig(k+1)];
	
	% Convergence test
	if max(abs(x(:,k+1)-x(:,k))) < eps
		break
	end
end