f = @(x,h) [1 + h^2 * (exp(x(2) .* sqrt(x(1))) + 3*x(1).^2); ...
            0.5 + h^2 * tan(exp(x(1)) + x(2).^2)];
J =  @(x, h) [(6*x(1)+0.5*exp(sqrt(x(1).*x(2))) .* sqrt(x(1))) (sqrt(x(1)).*exp(x(2).*sqrt(x(1)))); ...
              (exp(x(1)) .* sec(exp(x(1)) + x(2).^2).^2) (2*y .* sec(exp(x(1)) + x(2).^2).^2)];

x0 = [1;0.5];
for i = 1:100
    
end