%%%%%%%%%%%%%%%%%%%%%%%
%% A Pie in the Oven %%
%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
format long e


%% Setting
dt = 0.0001;
final_times = 1;
num_it=round(final_times/dt);
N = 10000; %Approximate number of nodes
Nbdr=300; %Number of nodes on the boundary
p = 3;% Polynomial Degree;
m = 7;% Order of the PHS --or-- shape parameter of GA
n = (p+2)*(p+1);% Number of neighbors (twice the number of polynomial terms)
P = (p+2)*(p+1)/2;% Number of polynomial terms
tol = 1e-12;

%Bounding box of the computational domain
R=1.5;
a1=-R; a2=R;
b1=-R; b2=R;

xfunc=@(t) ((17/31).*sin((235/57)+(-32).*t)+(19/17).*sin((192/55)+(-30).*t)+(47/32).*sin((69/25)+(-29).*t)+(35/26).*sin((75/34)+(-27).*t)+(6/31).*sin((23/10)+(-26).*t)+(35/43).*sin((10/33)+(-25).*t)+(126/43).*sin((421/158)+(-24).*t)+(143/57).*sin((35/22)+(-22).*t)+(106/27).*sin((84/29)+(-21).*t)+(88/25).*sin((23/27)+(-20).*t)+(74/27).*sin((53/22)+(-19).*t)+(44/53).*sin((117/25)+(-18).*t)+(126/25).*sin((88/49)+(-17).*t)+(79/11).*sin((43/26)+(-16).*t)+(43/12).*sin((41/17)+(-15).*t)+(47/27).*sin((244/81)+(-14).*t)+(8/5).*sin((79/19)+(-13).*t)+(373/46).*sin((109/38)+(-12).*t)+(1200/31).*sin((133/74)+(-11).*t)+(67/24).*sin((157/61)+(-10).*t)+(583/28).*sin((13/8)+(-8).*t)+(772/35).*sin((59/16)+(-7).*t)+(3705/46).*sin((117/50)+(-6).*t)+(862/13).*sin((19/8)+(-5).*t)+(6555/34).*sin((157/78)+(-3).*t)+(6949/13).*sin((83/27)+(-1).*t)+(-6805/54).*sin((1/145)+2.*t)+(-5207/37).*sin((49/74)+4.*t)+(-1811/58).*sin((55/43)+9.*t)+(-63/20).*sin((2/23)+23.*t)+(-266/177).*sin((13/18)+28.*t)+(-2/21).*sin((7/16)+31.*t))/1000;
yfunc=@(t) ((70/37).*sin((65/32)+(-32).*t)+(11/12).*sin((98/41)+(-31).*t)+(26/29).*sin((35/12)+(-30).*t)+(54/41).*sin((18/7)+(-29).*t)+(177/71).*sin((51/19)+(-27).*t)+(52/33).*sin((133/52)+(-26).*t)+(59/34).*sin((125/33)+(-26).*t)+(98/29).*sin((18/11)+(-25).*t)+(302/75).*sin((59/22)+(-24).*t)+(104/9).*sin((118/45)+(-22).*t)+(52/33).*sin((133/52)+(-21).*t)+(37/45).*sin((61/14)+(-20).*t)+(143/46).*sin((144/41)+(-19).*t)+(254/47).*sin((19/52)+(-18).*t)+(246/35).*sin((92/25)+(-17).*t)+(722/111).*sin((176/67)+(-16).*t)+(136/23).*sin((3/19)+(-15).*t)+(273/25).*sin((32/21)+(-13).*t)+(229/33).*sin((117/28)+(-12).*t)+(19/4).*sin((43/11)+(-11).*t)+(135/8).*sin((23/10)+(-10).*t)+(205/6).*sin((33/23)+(-8).*t)+(679/45).*sin((55/12)+(-7).*t)+(101/8).*sin((11/12)+(-6).*t)+(2760/59).*sin((40/11)+(-5).*t)+(1207/18).*sin((21/23)+(-4).*t)+(8566/27).*sin((39/28)+(-3).*t)+(12334/29).*sin((47/37)+(-2).*t)+(15410/39).*sin((185/41)+(-1).*t)+(-596/17).*sin((3/26)+9.*t)+(-247/28).*sin((25/21)+14.*t)+(-458/131).*sin((21/37)+23.*t)+(-41/36).*sin((7/8)+28.*t))/1000;
xpfunc =@(t) (-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t);
ypfunc =@(t) (-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t);

xfunc2=@(t) R*cos(t);
yfunc2=@(t) R*sin(t);
xpfunc2 =@(t) -R*sin(t);
ypfunc2 =@(t) R*cos(t);
%% End Settings %%%%%


%% RBF Definition --- PHS
phi =@(X1,X2,Y1,Y2,alpha) sqrt((X1-X2).^2+(Y1-Y2).^2).^alpha;
phi_x =@(X1,X2,Y1,Y2,alpha) -alpha*(X1-X2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_y =@(X1,X2,Y1,Y2,alpha) -alpha*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_xx =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(X1-X2).^2+(Y1-Y2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_yy =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(Y1-Y2).^2+(X1-X2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_lap =@(X1,X2,Y1,Y2,alpha) alpha^2*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);


%% Node Distribution
[x y] = hex(N,[a1,a2,b1,b2]);

%The nodes must lay inside the closed curve defined by the parametric
%equation {xfunc2(t(:)) yfunc2(t(:))}
%All other nodes are dismissed
t = linspace(0,2*pi,Nbdr+1); t = t(1:end-1);
dsites = [xfunc2(t(:)) yfunc2(t(:))];
scal = (xpfunc2(t(:)).^2+ypfunc2(t(:)).^2).^(1/2);
norm = [ypfunc2(t(:)) -xpfunc2(t(:))]./[scal scal];
[idx,dist] = knnsearch(dsites,[x y],'k',2);
node2curve=dsites(idx(:,2),:)-[x y];
node2curve_dot_norm=sum(norm(idx(:,2),:).*node2curve,2);
x(node2curve_dot_norm<=0)=[]; %Remove points outside the curve
y(node2curve_dot_norm<=0)=[];

x=[xfunc2(t(:));x];%Adds the points on curve as the boundary nodes
y=[yfunc2(t(:));y];
index_bdr = 1:Nbdr;

%% Initial condition
%%This is a step function (1 inside the curve and -1 outside)
t = linspace(0,2*pi,1000);
t = t(1:end-1);
dsites = [xfunc(t(:)) yfunc(t(:))];
scal = (xpfunc(t(:)).^2+ypfunc(t(:)).^2).^(1/2);
norm = [ypfunc(t(:)) -xpfunc(t(:))]./[scal scal];
%Decide whether a point is inside or outside the Pi
[idx,dist] = knnsearch(dsites,[x y],'k',2);
node2Pi=dsites(idx(:,2),:)-[x y];
node2Pi_dot_norm=sum(norm(idx(:,2),:).*node2Pi,2);
f = sign(node2Pi_dot_norm);


%% Nearest Neighbors
[idx,dist] = knnsearch([x y],[x y],'k',n);

%% Polynomial exponents
%Sets up help matrices - polynomial augmenting of radial functions
dsites = [x y];
x_exp = [];
y_exp = [];
for i_exp = 0:p
    x_exp = [x_exp 0:i_exp];
    y_exp = [y_exp i_exp-(0:i_exp)];
end

x_exp_Dx = x_exp-1;
x_exp_Dx(x_exp_Dx<0)=0;
x_exp_Dxx = x_exp-2;
x_exp_Dxx(x_exp_Dxx<0)=0;

y_exp_Dy = y_exp-1;
y_exp_Dy(y_exp_Dy<0)=0;
y_exp_Dyy = y_exp-2;
y_exp_Dyy(y_exp_Dyy<0)=0;



%% Initialize the differentiation matrices
Dx = sparse(N);
Dy = sparse(N);
Dlap = sparse(N);

%% Populating the Differentiation Matrices
for index = 1:length(x)
    
    xn = dsites(idx(index,:),1);
    yn = dsites(idx(index,:),2);
    
    [X1,X2]=meshgrid(xn);
    [Y1,Y2]=meshgrid(yn);
    
    A11 = phi(X1,X2,Y1,Y2,m);
    A12 = (xn(:,ones(1,P)).^x_exp(ones(1,n),:)).*(yn(:,ones(1,P)).^y_exp(ones(1,n),:));
	
    Bx11 = phi_x(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Bx12 = x_exp.*(xn(1,ones(1,P)).^x_exp_Dx).*(yn(1,ones(1,P)).^y_exp);
    
    By11 = phi_y(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    By12 = y_exp.*(xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp_Dy);
    
    Blap11 = phi_lap(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Blap12 = x_exp.*x_exp_Dx.*(xn(1,ones(1,P)).^x_exp_Dxx).*(yn(1,ones(1,P)).^y_exp)+y_exp.*y_exp_Dy.*(xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp_Dyy);
    
    A = [A11 A12;A12' zeros(P)];
    Bx = [Bx11 Bx12];
    By = [By11 By12];
    Blap = [Blap11 Blap12];
    
    Dx_local = Bx/A;
    Dy_local = By/A;
    Dlap_local = Blap/A;
    
    Dx(index,idx(index,:)) = Dx_local(1,1:n);
    Dy(index,idx(index,:)) = Dy_local(1,1:n);
    Dlap(index,idx(index,:)) = Dlap_local(1,1:n);
end

%Computes the time stepper - BDF1
D = speye(length(x))-dt*Dlap;
D(index_bdr,:)=0;
D(index_bdr,index_bdr)=speye(Nbdr); %Dirichlet BC


[X, Y] = meshgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100));

u=f;
for iter = 1:num_it
    
    u(index_bdr)=-1; %Dirichlet BC
    
    u= D\u; %Advances in time
    
	%Contour plot of solution
	figure(1)
	hold off
	[c h]=contourf(X,Y,griddata(x,y,u,X,Y),50);
	set(h, 'edgecolor','none');
	colormap(jet(256));
	colorbar;
	caxis([-1 1])
	axis([-R,R,-R,R])
	drawnow
	pause(.01)
	
	
	%{
%Plots the solution
    figure(1)
    hold off
    plot3(x,y,u,'.')
    title('Solving the Heat Equation using RBF-FD')
    view(-40, 38);
    axis equal; axis off;
    axis([-R, R, -R,R, -2,2])
    drawnow
	pause(.01)
	%}
end

% gmres(D,u,[],10^(-13),500);