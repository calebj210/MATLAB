function ier = halley(x0, f, df, ddf, maxIts, tol, root)
    fprintf("Halley's Method\n")
    for i = 1:maxIts
        F   = f(x0);
        DF  = df(x0);
        DDF = ddf(x0);
        
        denom = 2*DF.^2 - F*DDF;
        
        if denom == 0
            ier = 2;
            return
        end
        
        x1 = x0 - 2*F*DF / denom;
        
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