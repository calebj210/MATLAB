%% Function Definition for Backward Substitution
function x = BS(U,b)

[n,~]=size(b);
epsilon=10^(-10);

for j = n:-1:1
	if abs(U(j,j)) < epsilon           % Check if U is singular
		display("U is singular")
		return
	end
	
	b(j) = b(j)/U(j,j);                     % Solve for b_j
	for i = j-1:-1:1
		b(i) = b(i) - U(i,j)*b(j);         % Compute new b
	end
end
x=b;
end