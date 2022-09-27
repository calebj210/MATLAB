function [approx] = hypfun_F_trapz(a,b,c,z,N)
% function [approx] = hypfun_F_trapz(a,b,c,z,N)
% 
% Computes the hypergeometric function \mathrm{F}(a,b;c;z) 
% using an end-corrected trapezoidal quadrature.
%
% Author: Caleb Jacobs, Cecile Piret
% Date last modified: 26-Sep-2022

%% Beta function: Parameters
h  = 1 / (N + 1);
al = b - 1; 
be = c - b - 1;
f  = @(t) (1 - z.*t).^(-a);
ra = round(2 + 1 / (5*h));
rb = round(2 + 1 / (5*h));
n  = 4 * ra;

fa = @(Z) (exp(1i*pi*be)) * (f(Z) .* (1 - Z).^be);
fb = @(Z) f(Z) .* Z.^al;

t = linspace(0, 2*pi, n + 1)'; t(end)=[];
zsa = 0 + h*ra*exp(1i*t);
zsb = 1 + h*rb*exp(1i*t);

pa = ra;
pb = rb;

%Correction stencils
W_lft = stencil_circ(n, al,  1, pa, ra, exp(1i*t));
W_rt  = stencil_circ(n, be, -1, pb, rb, exp(1i*t));

zk = h*((pa) : (N - (pb - 1))); %TR path nodes

TR = h * sum(f(zk)    .* (zk.^al)   .* ((zk - 1).^be));
corr_a = sum(W_lft(:) .* fa(zsa(:))) * ( h).^(al + 1);
corr_b = sum(W_rt(:)  .* fb(zsb(:))) * (-h).^(be + 1);

approx = 1 / gamfun(b) / gamfun(c - b) * ...
        (TR + corr_a - corr_b);
end

function W = stencil_circ(n,gam,dir,p,ra,za)
%%Computes the endpoint correction stencils
% INPUTS
% n: Stencil Size (2n+1)x(2n+1)
% gam: Order of the fractrional derivative
% dir: direction of the path
% p: we include TR terms from k=1 to k=p-1
% za: stencil nodes
% OUTPUTS
% W: weights in the stencil
zk  = za(:).';      % Stencil nodes
lzk = length(zk);
ord = 0 : (lzk - 1);

%Calculate stencil weights
Rinv = diag(1 ./ (ra.^ord));
vS   = -(dir.^ord(:)) .* hurwitzZeta(-gam - ord(:), p);

W = 1 / lzk * (zk.^-ord(:)) * Rinv * vS;
end