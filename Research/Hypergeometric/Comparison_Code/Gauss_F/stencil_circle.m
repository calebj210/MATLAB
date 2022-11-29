function W = stencil_circle(al,dir,p,zk)
%%Computes the endpoint correction stencils
%INPUTS
%n: Stencil Size (2n+1)x(2n+1)
%al: Order of the fractrional derivative
%za: nodes in the stencil. Note: They are scaled by 1/h
%OUTPUTS
%W: weights in the stencil

%Calculate stencil weights
zk = zk(:).';
lzk=length(zk);
ord = (0:(lzk-1))';

A = zk.^(ord);
vS = (dir.^ord(:)).*(zeta(-al-ord(:))-sum(1./(1:(p-1)).^(-al-ord(:)),2));

cS = inv(A)*vS;

W = cS(:);

end