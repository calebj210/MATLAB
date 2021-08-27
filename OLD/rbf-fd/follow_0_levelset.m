clear all;
clc;
close all;

%---Inputs---
%xran = [-5 5]; yran = [-1 1]; x and y ranges
xran = [-1 1]; yran = [-1 1]; % x and y ranges
cN = 40;                      % number of points on curve
maxgN = 320;                 % max number of background nodes


%---Node Distribution---
box = [xran yran]; % bounding box

[x,y] = hex(cN,box); % hexagonal grid for background nodes
gN = length(x); % number of grid nodes

%Function and derivatives
%R = 5;
R = 1;
S =@(X,Y) sqrt(X.^2+Y.^2)-R^2;
Sx =@(X,Y) X./sqrt(X.^2 + Y.^2);
Sy =@(X,Y) Y./sqrt(X.^2 + Y.^2);

figure(1)
plot(x,y,'o')
axis equal
p = 10; %try p>1
disp('Step 1: Find nodes on the curve')

%%Charge (The highest values are associated to the densest regions)
%strength = 10; %for R = 5
strength = 1;   %for R = 1
Q = @(x,y) strength*ones(size(x));

%%Number of neighbours
n = 10;
delta = 10^(-4);
ux_old=ones(size(x));
uy_old=ones(size(x));

for count=1:100
    ux_old = x;
    uy_old = y;
    
    figure(1)
    hold off
    plot(x,y,'o')
    axis equal
    drawnow
    
    
    %%Looking for the n nearest neighbours
    zx = [x];
    zy = [y];
    
    [idx,dist]=knnsearch([zx zy],[x y],'k',n+1);
    
    rx_i = -zx(idx(:,2:end))+x(:,ones(1,n));
    ry_i = -zy(idx(:,2:end))+y(:,ones(1,n));

    Qi = Q(zx(idx(:,2:end)),zy(idx(:,2:end)));
   
    %Coulomb's law
    T = Qi.*Qi./(dist(:,2:end).^3);
    Fx = sum(T.*rx_i,2);
    Fy = sum(T.*ry_i,2);
    alpha = (S(x,y))./((Sx(x,y)).^2+(Sy(x,y)).^2);
    x = x+delta*Fx - alpha.*Sx(x,y);
    y = y+delta*Fy - alpha.*Sy(x,y);
end

gamma(:,1)=x;
gamma(:,2)=y;


% disp('Press a key')
% pause



disp('Step 2: Update nodes')

[x,y] = hex(maxgN,box); % hexagonal grid for background nodes

gN = length(x); % number of grid nodes

%%Charge (The highest values are associated to the densest regions)
%strength = 6; %for R = 5
strength = 1.2; %for R = 1
Q = @(x,y) strength*ones(size(x)); % exp(-10*((x-1/2).^2+(y-1/2).^2));

%%Number of neighbours
n = 10;
delta = 10^(-4);
temp=1;

while temp>6*10^(-3)
    
    figure(1)
    hold off
    plot(gamma(:,1),gamma(:,2),'ko')
    hold on
    plot(x,y,'*')
    axis equal
	xlim([-2,2])
	ylim([-2,2])
    drawnow
    
    %%Looking for the n nearest neighbours
    zx = [x;gamma(:,1)];
    zy = [y;gamma(:,2)];
    
    [idx,dist]=knnsearch([zx zy],[x y],'k',n+1);
    
    rx_i = -zx(idx(:,2:end))+x(:,ones(1,n));
    ry_i = -zy(idx(:,2:end))+y(:,ones(1,n));
    
    Qi = Q(zx(idx(:,2:end)),zy(idx(:,2:end)));
    
    %Coulomb's law
    T = Qi.*Qi./(dist(:,2:end).^3);
    Fx = sum(T.*rx_i,2);
    Fy = sum(T.*ry_i,2);
    alpha = 0.05*(S(x,y))./((Sx(x,y)).^2+(Sy(x,y)).^2);
    
    x_temp = x;
    y_temp = y; % new
    
    x = x+delta*Fx-alpha.*Sx(x,y);
    y = y+delta*Fy-alpha.*Sy(x,y);


    temp=[norm([x_temp;y_temp]-[x;y],inf)];
end