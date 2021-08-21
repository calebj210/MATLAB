%% Function Definition for Forward Substitution
function x = FS(L,b)

[n,~]=size(b);
epsilon=10^(-10);

for i = 1:n
	for j = 1:i-1
		b(i) = b(i) - L(i,j)*b(j);    % Compute new b
	end
	
	if abs(L(i,i)) < epsilon       % Check if L is singular
		display("L is singular")
		return
	end
	
	b(i) = b(i)./L(i,i);                % Solve for b_i
end
x=b;
end