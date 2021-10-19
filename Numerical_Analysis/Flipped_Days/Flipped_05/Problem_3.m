%% 
% Flipped day 5 problem
% Author: Caleb Jacobs
% Date last modified: 15-10-2021
clear
close

%% Intial data
f  = @(x) x ./ (1 + 25*x.^2);
df = @(x) (1 - 25*x.^2) ./ (1 + 25*x^2).^2;

N  = 21;

xi = linspace(-1, 1, N);
x  = linspace(-1, 1, 1000);

%% Standard polynomial using barycentric interpolation
w  = baryWeights(xi);
bf = BF(x, xi, w, f);

figure(1)
plot(x, f(x))
hold on
plot(x, bf)
title('Standard polynomial')

%% Hermite interpolant

%% Piecewise-linear continuous interpolant
lf = PL(x, xi, f);

figure(3)
plot(x, f(x))
hold on
plot(x, lf)
title('Piecewise-linear')

%% Piecewise-cubic Hermite interpolant
w  = pHWeights(xi, f, df);
hf = PH(x, xi, w);

figure(4)
plot(x, f(x))
hold on
plot(x, hf)
title('Piecewise-cubic Hermite')

%% Natural cubic spline

%% Functions
% Barycentric
function w = baryWeights(xi)
    n = length(xi);         % Number of nodes
    w = zeros(n, 1);        % Initialize weights
    
    % Compute weights
    for j = 1 : n
        lp = 1;
        for i = 1 : j - 1
            lp = lp * (xi(j) - xi(i));
        end
        for i = j + 1 : n
            lp = lp * (xi(j) - xi(i));
        end
        
        w(j) = 1 / lp;
    end
end

function func = BF(x, xi, w, f)
    n = length(x);          % Number of function evaluations
    N = length(xi);         % Number of nodes
    
    y = f(xi);              % Function value at each node
    
    func = zeros(n, 1);     % Function values
    for i = 1 : n
        num = 0;            % Numerator of barycentric form
        den = 0;            % Denominator of barycentric form
        for j = 1 : N
            if x(i) - xi(j) ~= 0
                tmp = w(j) / (x(i) - xi(j));
                num = num + tmp * y(j);
                den = den + tmp;
            else
                num = y(j);
                den = 1;
                break
            end
        end
        
        func(i) = num / den; % Compute function value
    end
end

% Hermite
function w = hWeights(x, f, df)
    w = 0;
end

% Piecewise linear
function func = PL(x, xi, f)
    n = length(x);              % Number of points to evaluate at
    N = length(xi);             % Number of nodes
    y = f(xi);                  % Different function evaluations
    
    func = zeros(n, 1);         % Initialize function values
    for i = 1:n
        for j = 2 : N
            if x(i) < xi(j)
                break
            end
        end
        
        func(i) = (y(j) - y(j - 1)) / (xi(j) - xi(j - 1)) * (x(i) - xi(j - 1)) + y(j - 1);
    end
end

% Piecewise Hermite 
function w = pHWeights(x, f, df)
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

function f = PH(x, xi, w)
    n = length(x);
    N = length(xi);
    f = zeros(n, 1);
    
    for j = 1 : n
        for i = 1 : N - 1
            if x(j) < xi(i + 1)
                break
            end
        end
        
        f(j) = [1 x(j) x(j).^2 x(j).^3] * w(:,i);
    end
end