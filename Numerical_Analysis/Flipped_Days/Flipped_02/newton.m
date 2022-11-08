function ier = newton(x0, f, df, maxIts, tol, root)
    fprintf("Newton's Method\n")
    for i = 1:maxIts
        denom = df(x0);
        
        if denom == 0
            ier = 2;
            return
        end
        
        x1 = x0 - f(x0) / denom;
        
        if abs(x1 - x0) <= tol
            ier = 0;
            fprintf('Root = %.16f\n\n', x1);
            return
        end
        
        x0 = x1;
        
        fprintf('%3d: %.16e\n', i, abs(x1 - root))
    end
    
    ier = 1;
    return
end