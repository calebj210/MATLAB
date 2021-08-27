% Time
time = 10;
iterations = 1000000;
t = linspace(0,time,iterations);
dt = time/iterations;

% Constants
Cd = .47;
roh = 1.2;
A = 10;
m = 5;
g = [0 -9.8];

% position vector
si = [0 0];
s = zeros(iterations,2);
s(1,:) = si;

% velocity vector
vi = [20 20];
v = zeros(iterations,2);
v(1,:) = vi;

% acceleration vector
a = zeros(iterations,2);

% Begin iterative simulation
for i = 1:iterations-1
	%{
	for j = 1:2 
		if v(i,j)<0
			a(i,j) = (((1/2)*roh*Cd*A)/m)*(v(i,j).^2)+g(j);
		else
			a(i,j) = -(((1/2)*roh*Cd*A)/m)*(v(i,j).^2)+g(j);
		end
	end
	%}
	vMag = sqrt((v(i,1)^2)+(v(i,2)^2));
	vUnit = (1/vMag)*v(i,:);
	a(i,:) = -(1/(2*m))*Cd*roh*(vMag^2)*A*vUnit+g;
	v(i+1,:) = v(i,:)+a(i,:)*dt;
	s(i+1,:) = s(i,:)+v(i,:)*dt;
end

subplot(3,2,1)
plot(s(:,1),s(:,2))
axis([0 inf 0 inf+.5])

subplot(3,2,3)
plot(t,v(:,1))
subplot(3,2,4)
plot(t,v(:,2))

subplot(3,2,5)
plot(t,a(:,1))
subplot(3,2,6)
plot(t,a(:,2))