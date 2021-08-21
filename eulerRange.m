% Base storage setup
xi = 0;
xf = 4;
yi = 2;
steps = 4;

dx = (xf-xi)/steps;
x = xi:dx:xf;
y = 1:size(x,2);

% Initial conditions
x(1) = xi;
y(1) = yi;

% Loop it boi
for i = 1:size(x,2)-1
	dy = x(i)-y(i);
	
	y(i+1) = y(i)+dy*dx;
end

plot(x,y)
grid on
