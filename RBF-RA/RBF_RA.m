%% Example 1: Interpolation problem using GA kernel
phi = @(e,r) exp(-(e*r).^2); % Gaussian
f = @(x,y) (1-(x.^2+y.^2)).*(sin(pi/2*(y-0.07))-0.5*cos(pi/2*(x+0.1)));
N = 60; hp = haltonset(2); y = 2*net(hp,N)-1; % Nodes/centers
M=2*N; x = 3/2*net(scramble(hp,'rr2'),M)-3/4; % Eval points
epsilon = linspace(0,1,101);

% Compute distances nodes/centers and eval points/centers
D = @(x,y)hypot(bsxfun(@minus,x(:,1),y(:,1)'),bsxfun(@minus,x(:,2),y(:,2)'));
ryy = D(y,y); rxy = D(x,y);

% Determine the radius
rad = fminbnd(@(e)norm(inv(phi(e,ryy)),inf)*norm(phi(1i*e,ryy),inf),0.1,20);
ryy = ryy*rad; % Re-scale the distances by the radius of the evaluation
rxy = rxy*rad; % contour so that a unit radius can be used.

% Compute the interpolant
rbfinterp=@(ep) phi(ep,rxy)*(phi(ep,ryy)\f(y(:,1),y(:,2)));
s = vvra(rbfinterp,epsilon/rad,1,64,64/4);

% Compute the difference between s and f and plot the results
error = max(abs(bsxfun(@minus,s,f(x(:,1),x(:,2)))));
semilogy(epsilon,error,'x-')
[minerr,pos] = min(error);
fprintf('Minimum error: %1.2e, Epsion: %1.2f\n',minerr,epsilon(pos));

%% Example 2: RBF-HFD weights for the 3-D Laplacian using IQ kernel
% Example for the standard 19 node stencil (6 implicit nodes) on lattice
% Explicit nodes
xhat = [[0,0,0];[-1,0,0];[1,0,0];[0,-1,0];[0,1,0];[0,0,-1];[0,0,1];...
[0,-1,-1];[0,-1,1];[0,1,-1];[0,1,1];[-1,0,-1];[-1,0,1];...
[1,0,-1];[1,0,1];[-1,-1,0];[-1,1,0];[1,-1,0];[1,1,0]];
N = size(xhat,1);

% Implicit nodes
yhat = [[-1,0,0];[1,0,0];[0,-1,0];[0,1,0];[0,0,-1];[0,0,1]];
plot3(xhat(:,1),xhat(:,2),xhat(:,3),'o',yhat(:,1),yhat(:,2),yhat(:,3),'.')

% The IQ kernel and it’s 3-D Laplacian and 3-D bi-harmonic
phi   = @(e,r) 1./(1 + (e*r).^2);
dphi  = @(e,r) 2*e^2*(-3 + (e*r).^2)./(1 + (e*r).^2).^3;
d2phi = @(e,r) 24*e^4*(5 + (e*r).^2.*(-10 + (e*r).^2))./(1 + (e*r).^2).^5;

% Compute the radius
rad = fzero(@(e)log10(cond(rbfhfd(e,xhat,yhat,phi,dphi,d2phi,1)))-6,[0.05,1]);
r   = sqrt(D(xhat,xhat).^2 + bsxfun(@minus,xhat(:,3),xhat(:,3)').^2);
rad = min(rad,0.95/max(r(:)));

% Re-scale the stencil nodes so that a unit evaluation radius can be used.
xhat = rad*xhat; yhat = rad*yhat;

% Compute the weights at various epsilon
epsilon = linspace(0,rad,11);
w = vvra(@rbfhfd,epsilon/rad,1,64,64/4,xhat,yhat,phi,dphi,d2phi);

% Undo the effects of re-scaling from the weights
w(1:size(xhat,1),:) = w(1:size(xhat,1),:)*rad^2;

% Compare the flat limit weights (epsilon=0) to the standard weights
ws = [-8,2/3*ones(1,6),ones(1,12)/3,-ones(1,6)/6]'; % standard weights
fprintf('Relative two norm difference: %1.2e\n',norm(w(:,1)-ws)/norm(ws));

function w = rbfhfd(ep,x,xh,phi,dphi,d2phi,flag)
    %RBFHFD Computes the RBF-HFD weights for the 3-D Laplacian.
    %
    % w = rbfhfd(Epsilon,X,Xh,Phi,DPhi,D2Phi) computes the RBF-HFD weights
    % at the given value of Epsilon and at the explicit stencil nodes X and
    % implicit (Hermite) stencil nodes Xh. Phi is a function handle for
    % computing the kernel phi(ep,r) used for generating the weights (e.g.
    % the Gaussian or inverse quadratic). Dphi and D2phi are function
    % handles for computing the 3-D Laplacian and bi-harmonic of Phi,
    % respectively.
    %
    % A = rbfhfd(Epsilon,X,Xh,Phi,DPhi,D2Phi,flag) returns the matrix for
    % computing the weights if the flag is non-zero, otherwise it returns
    % just the weights as described above.
    if nargin == 6
        flag = 0;                       % Return on the weights
    end
    N  = size(x,1);
    x  = [x;xh];
    ov = ones(1,size(x,1));

    % Compute the pairwise distances
    r = sqrt((x(:,1)*ov-(x(:,1)*ov).').^2 + ...
             (x(:,2)*ov-(x(:,2)*ov).').^2 + ...
             (x(:,3)*ov-(x(:,3)*ov).').^2);

    temp = dphi(ep,r(1:N,N+1:end));     % Construct the weight matrix
    A = [[phi(ep,r(1:N,1:N)) temp];[temp.' d2phi(ep,r(N+1:end,N+1:end))]];

    % Determine what needs to be returned.
    if flag == 0
        w = A \ [dphi(ep,r(1:N,1));d2phi(ep,r(N+1:end,1))];
    else
        w = A;
    end
end

function [R,b] = vvra(myfun,epsilon,rad,K,n,varargin)
    %VVRA Vector-valued rational approximation (VVRA) valid near the
    % origin to a vector-valued analytic function.
    %
    % [R,b] = vvra(myfun,Epsilon,Rad,K,n) generates a vector-valued rational
    % approximation to the vector-valued analytic function represented by
    % myfun and evaluates it at the values in Epsilon. The inputs are
    % described as follows:
    %
    % myfun is a function handle that, given a value of epsilon returns a
    % column vector corresponding to each component of the function
    % evaluated at epsilon.
    %
    % Epsilon is an array of values such that |Epsilon|<=Rad where the
    % rational approximation is to be evaluated.
    %
    % Rad is a scalar representing the radius of the circle centered at
    % the origin where it is numerically safe to evaluate myfun.
    %
    % K is the number of points to evaluate myfun at on the contour for
    % constructing the rational approximation. This number should be even.
    %
    % 2n is the degree of the common denominator to use in the
    % approximation.
    %
    % The columns in the output R represent the approximation to the
    % components of myfun at each value in the array Epsilon. The optional
    % output b contains the coefficients of the common denominator of the
    % vector-valued rational approximation.
    %
    % [R,b] = vvra(myfun,Epsilon,Rad,K,n,T1,T2,...) is the same as above, but
    % passes the optional arguments T1, T2, etc. to myfun, i.e.,
    % feval(myfun,epsilon,T1,T2,...).
    K = K+mod(K,2);                 % Force K to be even
    ang = pi/2*linspace(0,1,K+1)'; ang = ang(2:2:K);
    ei = rad*exp(1i*ang);           % The evaluation points (all in first quadrant)
    m = K-n;                        % Fix the degree of the numerator based on K and n
    W = feval(myfun,ei(1),varargin{:});
    M = numel(W);
    fv = zeros(K/2,1,M);
    fv(1,1,:) = W;

    % Loop over the evaluation points for F
    for k = 2:K/2
        W = feval(myfun,ei(k),varargin{:});
        fv(k,1,:) = W;
    end

    fmax = max(abs(fv),[],3);       % Find largest magnitude component for each k
    e = ei.^2;                      % Calculate the E matrix
    E = e(:,ones(1,m)); E(:,1) = 1./fmax; E = cumprod(E,2); % Scaled E-matrix
    f = E(:,1:n+1);                 % Create the F-matrices and RHS
    F = f(:,:,ones(1,M)).*fv(:,ones(1,n+1),:);
    g = F(:,1,:);
    F = -F(:,2:n+1,:);
    ER = [real(E);imag(E)];         % Separate E and F and the RHS g into real parts on
    FR = [real(F);imag(F)];         % top, and then the imag parts below
    gr = [real(g);imag(g)];
    [Q,R] = qr(ER); QT = Q';        % Factorize ER into Q*R
    R = R(1:m,:);                   % Remove bottom block of zeros from R

    % Update all FR matrices
    for k = 1:M
        FR(:,:,k) = QT*FR(:,:,k);
        gr(:,1,k) = QT*gr(:,1,k);
    end
    
    FT = FR(1:m,:,:);               % Separate F-blocks and g-blocks to the systems for
    FB = FR(m+1:K,:,:);             % determining the numerator and denominator coeffs.
    gt = gr(1:m,:,:);
    gb = gr(m+1:K,:,:);

    % Reshape these to be 2-D matrices
    FT = permute(FT,[1,3,2]); FT = reshape(FT,M*m,n);
    FB = permute(FB,[1,3,2]); FB = reshape(FB,M*(K-m),n);
    gt = permute(gt,[1,3,2]); gt = reshape(gt,M*m,1);
    gb = permute(gb,[1,3,2]); gb = reshape(gb,M*(K-m),1);
    b = FB\gb;                      % Obtain the coefficients of the denominator
    v = gt-FT*b; V = reshape(v,m,M);
    a = (R\V);                      % Obtain the coefficients of the numerators
    
    % Evaluate the rational approximations
    R = zeros(M,length(epsilon));
    b = [1;b];
    denomval = polyval2(b,epsilon);
    for ii = 1:M
        R(ii,:) = (polyval2(a(:,ii),epsilon)./denomval);
    end
end

function y = polyval2(p,x)
    %POLYVAL2 Evaluates the even polynomial
    % Y = P(1) + P(2)*Xˆ2 + ... + P(N)*Xˆ(2(N-1)) + P(N+1)*Xˆ2N
    % If X is a matrix or vector, the polynomial is evaluated at all
    % points in X (this is unlike the polyval function of matlab)
    y = zeros(size(x)); x = x.^2;

    for j=length(p):-1:1
        y = x.*y + p(j);
    end
end