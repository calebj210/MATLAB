function [h] = hypfun_F_trapz(a,b,c,z,N)
% function [h] = hypfun_F_trapz(a,b,c,z,N)
% 
% Computes the hypergeometric function \mathrm{F}(a,b;c;z)
% using a Pochhammer integral approach.
%
% Author: Caleb Jacobs, Cecile Piret
% Date last modified: 17-Oct-2022

% Integrand parameters
al = b - 1;
be = c - b - 1;
f = @(Z) (1 - Z).^(-a);

% Compute F using end correct trapezoidal rule
h = 1 / gamfun(b) / gamfun(c - b) / ...
    (z .^ (c - 1))                * ...
    hyp_integrate(0, z, f, al, be);
end