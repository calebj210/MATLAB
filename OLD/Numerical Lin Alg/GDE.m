function x = GDE(a,b)

%% Define tolerance and size of matrix %%
epsilon=10.^(-10);
[n,n]=size(a);

%% eliminate each column of A %%
for k=1:n
	if abs(a(k,k))<epsilon
		display("Matrix is singular or has a zero on the diagonal")
		return;
	end
	
	% Eliminate above diagonal %
	for i=1:k-1
		a(i,k)=-a(i,k)./a(k,k);
		for j=k+1:n
			a(i,j)=a(i,j)+a(k,j).*a(i,k);
		end
		b(i,1)=b(i,1)+b(k,1).*a(i,k); % Update b
	end
	
	% Eliminate below diagonal % 
	for i=k+1:n
		a(i,k)=-a(i,k)./a(k,k);
		for j=k+1:n
			a(i,j)=a(i,j)+a(k,j).*a(i,k);
		end
		b(i,1)=b(i,1)+b(k,1).*a(i,k); % Update b
	end
end

%% Solve Dx=b for x %%
for k=1:n
	b(k,1)=b(k,1)./a(k,k);
end

b
end