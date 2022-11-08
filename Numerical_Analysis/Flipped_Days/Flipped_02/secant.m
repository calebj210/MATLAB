function ier = secant(x0, x1, f, maxIts, tol, root)
    fprintf("Secant Method\n")
    for i = 1:maxIts
        denom = f(x1) - f(x0);
        
        if denom == 0
            ier = 2;
            return
        end
        
        x2 = x1 - f(x1) * (x1 - x0) / denom;
        
        if abs(x1 - x0) <= tol
            ier = 0;
            fprintf('Root = %.16f\n\n', x1);
            return
        end
        
        x0 = x1;
        x1 = x2;
        
        fprintf('%3d: %.16e\n', i, abs(x2 - root))
    end
    
    ier = 1;
    return
end