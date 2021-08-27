%Evenly spaced node distribution on unit circle.

close all
clear all

%---Inputs---
cN = 10000;
%---Parameterized Unit Circle---
xfunc = @(th) cos(th);
yfunc = @(th) sin(th);

t = linspace(0,2*pi,cN);

%storage for nodes on curve

gamma = zeros(cN,2);

gamma(:,1) = xfunc(t);
gamma(:,2) = yfunc(t);

figure
axis equal
hold on
plot(gamma(:,1),gamma(:,2),'bo')

