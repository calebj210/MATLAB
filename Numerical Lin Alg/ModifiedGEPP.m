%% Defining function to compute GE with PP given A and b %%
function [L,U,p,x,r,rnorm] = GEPP(a,b)

%% Setup %%
epsilon=10.^(-10);
[n,n]=size(a);
p=1:n;
singular=0;
A=a;
b0=b;

%% LU with PP %%
for k=1:n-1
	
	% Find max value in column under the diagonal and
	% permute the permutation vector
	[m,idx]=max(abs(a(p(k:end),k)));
	m=a(p(k+idx-1),k);
	temp=p(k+idx-1);
	p(k+idx-1)=p(k);
	p(k)=temp;
	
	% Flag if A is singular
	if abs(a(p(k),k))<epsilon
		display('Matrix A is singular')
		singular=1;
		continue;
	end
	
	% Compute parts of L and transform the required elements in A
	for i=k+1:n
		a(p(i),k)=a(p(i),k)./m;
		
		for j=k+1:n
			a(p(i),j)=a(p(i),j)-a(p(k),j).*a(p(i),k);
		end
	end
end

%% Check last entry if U is equal to 0 to see if A is singular
if abs(a(p(k+1),k+1))<epsilon
	display('Matrix A is singular')
	singular=1;
end

%% Column Oriented FS Solving Ly=b %%
if singular == 0
	for j=1:n
		for i=j+1:n
			b(p(i),1)=b(p(i),1)-a(p(i),j).*b(p(j),1);
		end
	end
end

%% Column Oriented BS Solving Ux=y %%
if singular == 0
	for j=n:-1:1
		b(p(j),1)=b(p(j),1)./a(p(j),j);
		for i=j-1:-1:1
			b(p(i),1)=b(p(i),1)-a(p(i),j).*b(p(j),1);
		end
	end
end

%% Formatting L,U, p, x, r, and ||r||_2 %%
% L %
L=eye(n);
for j=1:n-1
	for i=j+1:n
		L(i,j)=a(p(i),j);
	end
end

% U %
U=a(p(:),:)-L+eye(n);

% Solution vector x %
if singular == 0
	x=b(p(:),1);
else
	x=NaN;
end

% Permuation vector%
p=p';

% Residual and norm of residual %
if singular == 0
	r=b0-A*x;
	rnorm=norm(r);
else
	r=NaN;
	rnorm=NaN;
end

%% Display Results %%
L
U
p
x
r
rnorm
end