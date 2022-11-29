%Driver to test our method to evaluate hypergeometric functions
%Linked to my CP221006 handout.
close all
clear all

% Test: Compute 1F1(c;d;z)
c = 3.2;
d = 5.4;
z = 1+1i;

true = hypergeom(c, d, z)
approx = gamma(d)/(z^(d-1)*gamma(c)*gamma(d-c))*hyp_integrate(0,z,@(Z) exp(Z),c-1,d-c-1)

rel_error = norm(true-approx) / norm(true)