% Parameters
u0 = 1;
w  = 2;
v0 = 2;
x0 = [0; 0];

% Velocity field definition
V = @(x, t) [u0*sin(w * (t - x(2) / v0)) ; v0];

% Time constraints
dt = 0.01;
tf = 5;
T = 0:dt:tf;
xP1 = zeros(2, length(T));
xP50 = zeros(2, length(T));
xP200 = zeros(2, length(T));

% Initialize tracers
x = zeros(2, length(T));
for i = 1:length(T)
   x(:, i) = x0;
end

% Run simulation using Euler's Method
for i = 1:length(T)
    for j = 1:i
        xj = x(:, j);
        vdt = dt * V(xj, T(i));
        x(:, j) = xj + vdt;
    end
    
    xP1(:, i)   = x(:, 1);
    xP50(:, i)  = x(:, 50);
    xP200(:, i) = x(:, 200);
    
    % Plot streakline
    if mod(i, 10) == 0
        scatter(x(1,1:i), x(2,1:i),5,T(i:-1:1))
        hold on
        scatter(xP1(1,1:i), xP1(2,1:i), 5)
        scatter(xP50(1,1:i), xP50(2,1:i), 5)
        scatter(xP200(1,1:i), xP200(2,1:i), 5)
        hold off
        title('Streakline from Origin')
        axis equal
        xlim([-5,5])
        ylim([0,tf*v0])
        drawnow
    end
end