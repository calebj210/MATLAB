function f = rosenbrock(x)
    N = length(x);  % Number of dimensions
    
    f = 0;
    for i = 1:(N - 1)
        f = f + (100 * (x(i)^2 - x(i + 1))^2 + (x(i) - 1)^2);
    end
end

