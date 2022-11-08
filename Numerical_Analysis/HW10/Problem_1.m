%%
% Comparison of lagrange interpolant and boolean sum lagrange
% Author: Caleb Jacobs
% Date last modified: 05-11-2021

%% Settings
f = @(x, y) (x + y) ./ (1 + 25 * (x.^2 + y.^2));
% f = @(x, y) exp(x + y);

errs = zeros(6, 4);
N = [5 : 5 : 30];

%% Part a
for i = 1 : length(N)
n = N(i);
x = linspace(-1, 1, n);
y = linspace(-1, 1, n);
[xtmp, ytmp] = meshgrid(x, y);
z = f(xtmp, ytmp);

xx = linspace(-1, 1, 100);
yy = linspace(-1, 1, 100);
[X, Y] = meshgrid(xx, yy);
Z = f(X, Y);
ZI = evalLagrange(X, Y, x, y, z);
ZB = evalBool(X, Y, x, y, f); 

errs(i, 1:2) = [norm(abs(Z(:) - ZI(:)), 2), norm(abs(Z(:) - ZB(:)), 2)];

%% Generate error plots
f1 = figure(1);
surf(X, Y, abs(Z - ZI))
title(sprintf('Equispaced Lagrange Absolute Error n = %d', n))

f2 = figure(2);
surf(X, Y, abs(Z - ZB))
title(sprintf('Equispaced Boolean Lagrange Absolute Error n = %d', n))

%% Part b
% n = 30;
k = 1 : n;
x = cos((2*k - 1) * pi / (2 * n));
y = cos((2*k - 1) * pi / (2 * n));
[xtmp, ytmp] = meshgrid(x, y);
z = f(xtmp, ytmp);

xx = linspace(-1, 1, 100);
yy = linspace(-1, 1, 100);
[X, Y] = meshgrid(xx, yy);
Z = f(X, Y);
ZI = evalLagrange(X, Y, x, y, z);
ZB = evalBool(X, Y, x, y, f);

errs(i, 3:4) = [norm(abs(Z(:) - ZI(:)), 2), norm(abs(Z(:) - ZB(:)), 2)];

f3 = figure(3);
surf(X, Y, abs(Z - ZI))
title(sprintf('Chebyshev Lagrange Absolute Error n = %d', n))

f4 = figure(4);
surf(X, Y, abs(Z - ZB))
title(sprintf('Chebyshev Boolean Lagrange Absolute Error n = %d', n))

%% Save figures
pause

% saveas(f1, sprintf('images/EL%02d.png', n))
% saveas(f2, sprintf('images/CL%02d.png', n))
% saveas(f3, sprintf('images/EB%02d.png', n))
% saveas(f4, sprintf('images/CB%02d.png', n))
end

figure(10)
semilogy(N, errs, 'LineWidth', 3)
legend('Equispaced Lagrange', ...
       'Equispaced Boolean', ...
       'Chebyshev Lagrange', ...
       'Chebyshev Boolean')
title('Absolute Error of Different Interpolants')
xlabel('Number of nodes')
ylabel('Absolute error')


%% Lagrange basis
function val = L(x, xi, j)
    n = length(xi);
    
    val = 1;
    for i = 1 : j - 1
        val = val * (x - xi(i)) / (xi(j) - xi(i));
    end
    for i = j + 1 : n
        val = val * (x - xi(i)) / (xi(j) - xi(i));
    end
end

%% Evaluate lagrange interpolation
function Z = evalLagrange(X, Y, x, y, z)
    Z = zeros(size(X));
    X = X(:);
    Y = Y(:);
    
    for i = 1 : length(X)
        for j = 1 : length(y)
            for k = 1 : length(x)
               Z(i) = Z(i) + z(j, k) * L(X(i), x, k) * L(Y(i), y, j);
            end
        end
    end
end

%% Evaluate boolean lagrange interpolation
function Z = evalBool(X, Y, x, y, f)
    p = length(x);
    q = length(y);
    Z = zeros(size(X));
    X = X(:);
    Y = Y(:);
    
    for i = 1 : length(X)
        sx = 0;
        for j = 1 : p
            sx = sx + f(x(j), Y(i)) * L(X(i), x, j);
        end
        
        sy = 0;
        for k = 1 : q
            sy = sy + f(X(i), y(k)) * L(Y(i), y, k);
        end
        
        cor = 0;
        for j = 1 : p
            for k = 1 : q
                cor = f(x(j), y(k)) * L(X(i), x, j) * L(Y(i), y, k);
            end
        end
        
        Z(i) = sx + sy - cor;
    end
end