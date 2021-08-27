%% Code for 6.a %%

% L=[1,0,0; 4, 10.^(-5),0; 6,7,-1]; % Define L for 5.a
% b=[1,-2,0];                              % Define b for 5.a
[c n]=size(b);
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
y1 = b


%% Code for 6.b %%
L=[-10,0,0,0; 1,1,0,0; -2,1,-6,0; 1,1,-1,3]; % Define L for 5.b
b=[0,0,-1,2];                                          % Define b for 5.b
[c n]=size(b);
epsilon=10^(-10);

for i = 1:n
	for j = 1:i-1
		b(i) = b(i) - L(i,j)*b(j);                      % Compute new b
	end
	
	if abs(L(i,i)) < epsilon                        % Check if L is singular
		display("L is singular")
		return
	end
	
	b(i) = b(i)/L(i,i);                                % Solve for b_i
end
y2 = b