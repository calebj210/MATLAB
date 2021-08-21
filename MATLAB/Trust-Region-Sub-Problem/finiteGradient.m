% Finite differenced gradient
function g = finiteGradient(f, x)
    D = length(x);       % Number of dimensions
    h = nthroot(eps, 3); % Optimal step length
    
    g = zeros(D,1);      % Initialize the gradient vector
    
    % Populate the gradient vector
    for i = 1:D
        forward  = x;           % Forward step
        backward = x;           % Backward step
        forward(i)  = x(i) + h;
        backward(i) = x(i) - h;
        
        g(i) = (f(forward) - f(backward)) / (2*h);
    end
end