%% Use LSSolve to find LS solution %%
close all
clear
format shortE
format loose

% Store given data
y = [0.0; 1.0; 2.0; 3.0; 4.0; 5.0];
b = [2.5; 9.2; 13.3; 26.7; 31.8; 50.4];

% Define A matrix for quadratic model
A = zeros(6,3);
for j = 1:3 
	A(:,j) = y.^(j-1);
end

% Display Initial A and b
disp('**********************************************')
disp('Initial A and b matrices')
disp('**********************************************')
A
b

% Solve LS probelm Ax =~ b
x = LSSolve(A,b);

% List coefficients individually to quadtratic model
disp('Where the coefficents to our quadratic model, f(y) = x1+x2*y+x3*y^2 are: ')
disp(' ')
disp(['x1 = ' num2str(x(1),15)])
disp(['x2 = ' num2str(x(2),15)])
disp(['x3 = ' num2str(x(3),15)])

% Plot data and solution
scatter(y,b,'filled')
hold on

t = linspace(0,5,100);
tt = x(1) + x(2).*t + x(3).*t.^2;
plot(t,tt,'LineWidth', 3)

title('Data and Solution to LS Problem')
legend('Data','x1+x2*y+x3*y^2')
xlabel('y')
ylabel('T(y)')
grid on