function f = rastrigin(x)
    D = length(x);
    f = 10 * D;
    for i = 1:D
        f = f + x(i).^2 - 10 * cos(2 * pi * x(i));
    end
end

