%% Function Definition for LSSolve, Least Squares Solver via QR with HH %%

function x = LSSolve(a,b)

% Settings
epsilon = 10.^(-10);

% Find dimensions of A
[m,n] = size(a);

%% Reduce A and b using QR %%
for k = 1:n
	alpha = -sign(a(k,k)).* norm(a(k:m,k));
	vk = [zeros(k-1,1); a(k,k) - alpha; a(k+1:m,k)];
	
    % Compute vk'*vk;
	beta = 0;
	for i = k:m
		beta = beta + vk(i).^2;
	end
	
	% Continue if beta is nearly zero to avoid OF
	if beta < epsilon
		continue
	end
	
    % Update A
	a(k,k) = alpha;
	a(k+1:m,k) = 0;
	
	for j = k+1:n
        % Compute vk'*a_ j;
		gamma = 0;
		for i = k:m
			gamma = gamma + vk(i).*a(i,j);
		end
		
		s = 2.*gamma./beta;
		
		for i = k:m
			a(i,j) = a(i,j) - s.*vk(i);
		end
	end
	
	%% Update b %%
    % Compute vk'*b;
	gamma = 0;
		for i = k:m
			gamma = gamma + vk(i).*b(i);
		end
	s = 2.*gamma./beta;
	
	for i = k:m
		b(i) = b(i) - s.*vk(i);
	end
	
	% Results of factorization
	disp('**********************************************')
	disp(['vk and reduced A and b for step  ', num2str(k)])
	if k == n
		disp('which is also the final reduced A and b:')
	end
	disp('**********************************************')
	
	vk
	a
	b
end

%% Solve Rx = c1 via BS %%
x = BS(a(1:n,:),b(1:n));

% Display solution x
disp('***********')
disp('Solution')
disp('***********')
x

% Compute ||r||^2
r = 0;
for i = n+1:m
	r = r + b(i).^2;
end

disp('||r||^2 =')
disp(' ')
disp(r)
end









