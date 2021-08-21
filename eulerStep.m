% Base storage setup
dx = 1;
x = 0:dx:10000;
y = 1:size(x,2);

% Initial conditions
x(1) = 0;
y(1) = 1000;

% Loop it boi
for i = 1:size(x,2)-1
	dy = .0017*y(i)*(1-(y(i)/5900));
	
	y(i+1) = y(i)+dy*dx;
end

plot(x,y)
