function W = stencil_circ(n,gam,dir,p,zk)
%%Computes the endpoint correction stencils
%INPUTS
%n: Number of nodes in the stencil 
%gam: Exponent of the singular term at the integral end
%dir: direction of the path
%p: p=1,2,3,... Shrinks the trapezoidal rule segment so that it only starts
%p nodes away from the integration path end.
%zk: stencil nodes -- we assume here that they belong to a circle with n
%equispaced nodes, starting with zk(1)=R, where R is the radius of the
%circle
%OUTPUTS
%W: weights in the stencil
zk = zk(:).'; % Stencil nodes
lzk=length(zk);
ord = 0:(lzk-1);
R = zk(1);% R is the radius of the correction stencil

%Calculate stencil weights
Rinv = diag(1./(R.^ord));
vS = hurwitzZeta(-gam-ord(:),p).*(dir.^ord(:));
RHS = Rinv*vS;
W = 1/lzk*fft(RHS);
end