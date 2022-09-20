function [h] = hypfun_F_pochham(a,b,c,Z,N)
% function [h] = hypfun_F_pochham(a,b,c,z,N)
% 
% Computes the hypergeometric function \mathrm{F}(a,b;c;z)
% using a Pochhammer integral approach.
%
% Author: Caleb Jacobs, Cecile Piret
% Date last modified: 19-Sep-2022


% Define base integrand
al = b - 1;
be = c - b - 1;
f = @(t) (t.^al) .* (1 - t).^be .* (1 - Z*t).^(-a);
x = @(t) -4/5*sin(2*t) + 1/2;
y = @(t) 1/10*sin(t) - 1/5*sin(3*t) - 1/10*sin(5*t);
z = @(t) x(t) + 1i*y(t);

% Computing dz for the trapezoidal rule
dx = @(t) -8/5*cos(2*t);
dy = @(t) 1/10*cos(t) - 3/5*cos(3*t) - 1/2*cos(5*t);
dz = @(t) dx(t) + 1i*dy(t);

% Computing the function values on the correct sheets
t  = pi/4+linspace(0,2*pi,N+1)'; t(1)=[];
t1 = t(t>pi/4 & t<=3*pi/4);
t2 = t(t>3*pi/4 & t<=5*pi/4);
t3 = t(t>5*pi/4 & t<=7*pi/4);
t4 = t(t>7*pi/4 & t<=9*pi/4+eps);
x1 = x(t1); x2 = x(t2); x3 = x(t3); x4 = x(t4);
y1 = y(t1); y2 = y(t2); y3 = y(t3); y4 = y(t4);
z1 = z(t1); z2 = z(t2); z3 = z(t3); z4 = z(t4); 
f1 = f(z1);
f2 = exp(-1i*2*pi*be)*f(z2);
f3 = exp(1i*2*pi*al)*exp(-1i*2*pi*be)*f(z3);
f4 = exp(1i*2*pi*al)*f(z4);

% Compute Pochhammer contour integral
h = 1 / gamfun(b) / gamfun(c - b) * ...
    1 / (1 - exp( 1i*2*pi*al))    * ...
    1 / (1 - exp(-1i*2*pi*be))    * ...
    (2*pi/N) * sum([f1(:); f2(:); f3(:); f4(:)] .* dz(t));
end