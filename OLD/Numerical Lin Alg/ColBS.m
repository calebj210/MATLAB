%% Code for 5.a %%

U=[-10,1,-1,0; 0, 10.^(-2),1,-1; 0,0,10.^2,1; 0,0,0,10.^(-3)]; % Define U for 4.a
b=[1,2,3,4];                                                                   % Define b for 4.a
[c n]=size(b);
epsilon=10^(-10);

for j = n:-1:1
	if abs(U(j,j)) < epsilon                                               % Check if U is singular
		display("U is singular")
		return
	end
	
	b(j) = b(j)/U(j,j);                                                          % Solve for b_j
	for i = j-1:-1:1
		b(i) = b(i) - U(i,j)*b(j);                                              % Compute new b
	end
end
x1=b


%% Code for 5.b %%

U=[2,-1,8;0,.1,4; 0,0,-1];             % Define U for 4.b
b=[1,0,3];                                  % Define b for 4.b
[c n]=size(b);
epsilon=10^(-10);

for j = n:-1:1
	if abs(U(j,j)) < epsilon           % Check if U is singular
		display("U is singular")
		return
	end
	
	b(j) = b(j)/U(j,j);                    % Solve for b_j
	for i = j-1:-1:1
		b(i) = b(i) - U(i,j)*b(j);        % Compute new b
	end
end
x2=b