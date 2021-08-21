%% Defining function to compute LU factorization with PP given A%%
function [L,U,P] = LU(a)

%% Setup %%
epsilon=10.^(-10);
[n,n]=size(a);
p=1:n;

%% LU with PP %%
for k=1:n-1
	
	% Find max value in column under the diagonal and
	% permute the permutation vector
	[m,idx]=max(abs(a(p(k:end),k)));
	m=a(p(k+idx-1),k);
	temp=p(k+idx-1);
	p(k+idx-1)=p(k);
	p(k)=temp;
	
	% Skip if A is singular
	if abs(a(p(k),k))<epsilon
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


% Permuation Matrix %
P=eye(n);
P=P(p,:);

%% Display Results %%
L
U
P
end