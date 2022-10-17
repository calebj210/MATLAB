%Driver to test our method to evaluate hypergeometric functions
%Linked to my CP221006 handout.
close all
clear all

h=0.4;
nxy = [-2,2,-2,2];
x = (nxy(1):h:nxy(2));
y = (nxy(3):h:nxy(4));
[xr,xi] = meshgrid(x,y(end:-1:1)); z = xr+1i*xi;

% Test1: Compute 1F1(c;d;z)
a = .2;
c = -2.2;
d = -4.4;

for i=1:length(z(:))
    true(i) = hypergeom([a c], d, z(i));
    approx(i) = gamma(d)/(z(i)^(d-1)*gamma(c)*gamma(d-c))*hyp_integrate(0,z(i),@(Z) 1./(1-Z).^a,c-1,d-c-1);
    rel_error(i) = norm(true-approx)/norm(true);
end

subplot(3,2,1)
mesh(xr,xi,reshape(real(true),size(xr)))
title(['Re _2F_1(' num2str(a) ',' num2str(c) ';' num2str(d) ';z)'])

subplot(3,2,2)
mesh(xr,xi,reshape(imag(true),size(xr)))
title(['Im _2F_1(' num2str(a) ',' num2str(c) ';' num2str(d) ';z)'])

subplot(3,2,3)
mesh(xr,xi,reshape(real(approx),size(xr)))
title('Approximation of the real part')

subplot(3,2,4)
mesh(xr,xi,reshape(imag(approx),size(xr)))
title('Approximation of the imaginary part')

subplot(3,1,3)
mesh(xr,xi,reshape(log10(abs(rel_error)),size(xr)))
title('log10(Relative Error)')