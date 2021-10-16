%% 
% Flipped day 5 problem
% Author: Caleb Jacobs
% Date last modified: 15-10-2021
f  = @(x) x ./ (1 + 25*x.^2);
df = @(x) (1 - 25*x.^2) ./ (1 + 25*x^2).^2;

N  = ;

xi = linspace(-1, 1, N);
x  = linspace(-1, 1, 100);
w  = hermiteWeights(xi, f, df);
hf = zeros(N, 1);
for i = 1 : length(x)
    hf(i) = H(x(i), xi, w);
end

figure(1)
plot(x, f(x))
hold on
plot(x, hf)

figure(2)
semilogy(x, abs(f(x)' - hf))

function w = hermiteWeights(x, f, df)
    n = length(x);
    w = zeros(4, length(x) - 1);
    
    p  = @(x) [1 x x.^2 x.^3];
    dp = @(x) [0 1 2*x 3*x.^2];
    
    for i = 1 : n - 1
        A = [p(x(i)); p(x(i+1)); dp(x(i)); dp(x(i+1))];
        b = [f(x(i)); f(x(i+1)); df(x(i)); df(x(i+1))];
        
        w(:,i) = A \ b;
    end
end

function f = H(x, xi, w)
    for i = 1 : length(xi) - 1
        if x < xi(i + 1)
            break
        end
    end
    
    f = [1 x x.^2 x.^3] * w(:,i);
end