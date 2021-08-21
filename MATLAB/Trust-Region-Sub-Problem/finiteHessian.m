function H = finiteHessian(f, x)
    D = length(x);          % Number of dimensions

    H = diag(df_xx(f, x));  % Initialize diagonal elements of Hessian
    
    for i = 1:D
        for j = (i + 1):D
            H(i,j) = df_xy(f, x, i, j);
            H(j,i) = H(i,j);
        end
    end
end


function df = df_xy(f, x, i, j)
    D = length(x);       % Number of dimensions
    h = nthroot(eps, 4); % Optimal step size
    
    hi = zeros(D, 1);    % ith step vector
    hj = zeros(D, 1);    % jth step vector
    hi(i) = h;
    hj(j) = h;
    
    df = (f(x + hi + hj) - f(x + hi - hj) - f(x - hi + hj) + f(x - hi - hj)) / (4 * h^2);
end

function df = df_xx(f, x)
    D = length(x);       % Number of dimensions
    h = nthroot(eps, 4); % Optimal step size
    
    df = zeros(D,1);     % Initialize the second derivative vector
    
    % Populate the gradient vector
    for i = 1:D
        forward  = x;           % Forward step
        backward = x;           % Backward step
        forward(i)  = x(i) + h;
        backward(i) = x(i) - h;
        
        df(i) = (f(forward) - 2*f(x) + f(backward)) / (h^2);
    end
end