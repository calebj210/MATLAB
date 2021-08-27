function v = myDE(z,Dx,Dy)
    gradz = [Dx*z' Dy*z'];
    
    for i = 1: length(z)
        normals(i,:) = gradz(i,:)./norm(gradz(i,:));    
    end
    
    kappa = abs(Dx*normals(:,1) + Dy*normals(:,2));
    
    for i = 1:length(z)
        v(i,1) = kappa(i).*(normals(i,1)*gradz(i,1) + normals(i,2)*gradz(i,2)); 
    end
end