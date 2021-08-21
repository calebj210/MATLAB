%% Defining function to compute GE with PP given A and b %%
function [x,L,U,P] = GEPP(a,b)

%% Setup %%
epsilon=10.^(-10);
[~,n]=size(a);
p=1:n;
singular=0;

%% LU with PP %%
for k=1:n-1
	
	% Find max value in column under the diagonal and
	% permute the permutation vector
	[~,idx]=max(abs(a(p(k:end),k)));
	m=a(p(k+idx-1),k);
	temp=p(k+idx-1);
	p(k+idx-1)=p(k);
	p(k)=temp;
	
	% Flag if A is singular
	if abs(a(p(k),k))<epsilon
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
	singular=1;
end
if singular ~= 0
	disp('Matrix A is Singular')
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

%% Formatting L,U, P, x%%
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
P=eye(n);
P=eye(n)(p,:);
end