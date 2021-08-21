function B = getEllipticalMatrix(d, scale)
    N = length(d);          % Number of dimensions
    normd = d / norm(d);
    
    % Perform Givens Rotations until matrix is aligned with d
    B = gramschmidt(normd); % Get orthogonal basis containing d
    D = eye(N);             % Initialize eigenvalue matrix
    D(1,1) = scale;         % Set eigenvalue for d
    
    B = B * D * B';         % Compute elliptical matrix
end

% Simple gramschmidt routine
function [U] = gramschmidt(d)
    N = length(d);
    U = rand(N);
    U(:,1) = d;
    for i = 2:N
        for j = 1:i-1
            U(:,i) = U(:,i) - (U(:,j)' * U(:,i)) / (norm(U(:,j)))^2 * U(:,j);
        end
        U(:,i) = U(:,i) / norm(U(:,i));
    end
end