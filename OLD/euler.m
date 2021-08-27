% Base storage setup
dx = .0001;
x = 0:dx:1;
y = 1:size(x,2);

% Initial conditions
x(1) = 0;
y(1) = 6;

% Loop it boi
for i = 1:size(x,2)-1
	dy = (15*x(i)^2)-(3*(x(i)^2)*y(i));
	
	y(i+1) = y(i)+dy*dx;
end

plot(x,y)
