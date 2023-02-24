function [g] = gamfun(z)
% function [g] = gamfun(z)
% 
% Master code for computing the Gamma function Gamma(z).
% 
% Uses Godfrey (Lanczos) and Stirling series to compute Gamma function for
% different parameter regimes.
% 
% Reliable up to |z|=171, at which point numerical overflow occurs.
% 
% Copyright John W. Pearson 2014

g = zeros(size(z));
for i = 1 : length(z)
    if real(z(i))<0.5 % transformation
        g(i) = pi/(sin(pi*z(i))*gamfun(1-z(i)));
    elseif imag(z(i))==0 && real(z(i))==fix(real(z(i))) % integer case
        g(i) = factorial(z(i)-1); % use built-in MATLAB routine
    elseif abs(z(i))>80
        g(i) = gamfun_stirlingbernoulli(z(i),eps);
    else
        g(i) = gamfun_godfrey(z(i));
    end
end
