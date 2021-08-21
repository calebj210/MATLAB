%%% Level-Set Methodology - Merging Circles

close all
clear all
format long e

%% Setting
dt = 0.00001;
final_times = 1;
time_init = 0;
num_it=round(final_times/dt);
N = 1000; %Approximate number of background nodes
%Nbdr=300; %Number of nodes on the boundary
p = 3;% Polynomial Degree;
m = 7;% Order of the PHS --or-- shape parameter of GA
n = (p+2)*(p+1);% Number of neighbors (twice the number of polynomial terms)
P = (p+2)*(p+1)/2;% Number of polynomial terms
tol = 1e-12;

%% Bounding box of the computational domain
R=0.3;
a1=-R; a2=R;
b1=-R; b2=R;
radius = 0.13;
aa1 = 0.1;
bb1 = 0.1;

aa2 = -0.1;
bb2 = -0.1;

%% RBF Definition --- PHS
phi =@(X1,X2,Y1,Y2,alpha) sqrt((X1-X2).^2+(Y1-Y2).^2).^alpha;
phi_x =@(X1,X2,Y1,Y2,alpha) -alpha*(X1-X2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_y =@(X1,X2,Y1,Y2,alpha) -alpha*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_xx =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(X1-X2).^2+(Y1-Y2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_yy =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(Y1-Y2).^2+(X1-X2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_xy =@(X1,X2,Y1,Y2,alpha) alpha*(alpha-1)*(X1-X2).*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_lap =@(X1,X2,Y1,Y2,alpha) alpha^2*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);

        %% RBF PHS 1D Definitions ---------------------------------------------------------
       phi1D =@(X1,X2,m) (abs(X1-X2)).^m;
       phi1D_Dx = @(X1,X2,m) -m.*sign(X1-X2).*((abs(X1-X2)).^(m-1));
       phi1D_Dxx = @(X1,X2,m) m.*(m-1).*((abs(X1-X2)).^(m-2));

        %% Polynomial Exponents 1D
        x1D_exp = [];
        for i_exp = 0:p
            x1D_exp=[x1D_exp i_exp];
        end
        x1D_exp_Dx = x1D_exp-1;
        x1D_exp_Dx(x1D_exp_Dx<0)=0;
        x1D_exp_Dxx = x1D_exp-2;
        x1D_exp_Dxx(x1D_exp_Dxx<0)=0;
        

%% Node Distribution
[x y] = hex(N,[a1,a2,b1,b2]);

%% Initial condition
%%This is a step function (1 inside the curve and -1 outside)

dsites = [x y];

for index = 1:length(dsites)
if dsites(index,1) >= -dsites(index,2)
init_distance(index) = sqrt((dsites(index,1) - aa1).^2+(dsites(index,2) - bb1).^2) - radius;
elseif dsites(index,1) < -dsites(index,2)
  init_distance(index) = sqrt((dsites(index,1) - aa2).^2+(dsites(index,2) - bb2).^2) - radius;
end
end
%initial zero level set (zls)
zlsfuncx1 =@(t) radius.*cos(t) + aa1;
zlsfuncy1 =@(t) radius.*sin(t) + bb1;

zlsfuncx2 =@(t) radius.*cos(t) + aa2;
zlsfuncy2 =@(t) radius.*sin(t) + bb2;

z = init_distance;


%% Initialize zero level set
% Number of Nodes
cN = 512;                      % number of points on both circles
cN = cN +2;                   % eliminate repeat end node
gamma = zeros(cN,2);

t = linspace(0,2*pi,cN/2)';

gamma(1:(cN/2),1) = zlsfuncx1(t);
gamma(1:(cN/2),2) = zlsfuncy1(t);

gamma((cN/2+1):end,1) = zlsfuncx2(t);
gamma((cN/2+1):end,2) = zlsfuncy2(t);


gamma(cN/2+2,:) = [];
gamma(cN/2,:) = [];


% % figure
% % % Plot initial zls and background nodes
% % hold on
% % axis equal
% % plot(x,y,'bo')
% % plot(gamma(:,1),gamma(:,2),'rx')

%% Nearest Neighbors
[idx_band,dist_band] = knnsearch([x y],[x y],'k',n); %KNN amongst band nodes


%% Polynomial exponents
%Sets up help matrices - polynomial augmenting of radial functions
x_exp = [];
y_exp = [];
for i_exp = 0:p;
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

exp_Dxy = x_exp.*y_exp;
x_exp_Dxy = exp_Dxy-y_exp_Dy;
x_exp_Dxy = x_exp_Dxy - 1;
x_exp_Dxy(x_exp_Dxy<0) = 0;
y_exp_Dxy = exp_Dxy-x_exp_Dx;
y_exp_Dxy = y_exp_Dxy -1;
y_exp_Dxy(y_exp_Dxy<0) = 0;


%% Initialize the differentiation matrices
Dx = sparse(N);
Dy = sparse(N);
Dlap = sparse(N);

for time = 0:num_it


    
%% Populating the Differentiation Matrices
for index = 1:length(x)
    
    
    xn = dsites(idx_band(index,:),1);
    
    yn = dsites(idx_band(index,:),2);
    
%%% What does this piece of code do?
    xn = xn - xn(1,1).*ones(length(xn),1);
    yn = yn - yn(1,1).*ones(length(yn),1);
    
    [X1,X2]=meshgrid(xn);
    [Y1,Y2]=meshgrid(yn);
    
    A11 = phi(X1,X2,Y1,Y2,m);
    A12 = (xn(:,ones(1,P)).^x_exp(ones(1,n),:)).*(yn(:,ones(1,P)).^y_exp(ones(1,n),:));
    
    Ax11 = phi_x(X1,X2,Y1,Y2,m);
    Ax12 = x_exp.*(xn(:,ones(1,P)).^x_exp_Dx).*(yn(:,ones(1,P)).^y_exp);
    
    Axx11 = phi_xx(X1,X2,Y1,Y2,m);
    Axx12 = x_exp.*x_exp_Dx.*(xn(:,ones(1,P)).^x_exp_Dxx).*(yn(:,ones(1,P)).^y_exp);
    
    Ay11 = phi_y(X1,X2,Y1,Y2,m);
    Ay12 = y_exp.*(xn(:,ones(1,P)).^x_exp).*(yn(:,ones(1,P)).^y_exp_Dy);
    
    Ayy11 = phi_yy(X1,X2,Y1,Y2,m);
    Ayy12 = y_exp.*y_exp_Dy.*(xn(:,ones(1,P)).^x_exp).*(yn(:,ones(1,P)).^y_exp_Dyy);
    
    Axy11 = phi_xy(X1,X2,Y1,Y2,m);
    Axy12 = x_exp.*y_exp.*(xn(:,ones(1,P)).^x_exp_Dxy).*(yn(:,ones(1,P)).^y_exp_Dxy);
    
    
    B11 = z(idx_band(index,:));
    B12 = (xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp);
    
    Bx11 = phi_x(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Bx12 = x_exp.*(xn(1,ones(1,P)).^x_exp_Dx).*(yn(1,ones(1,P)).^y_exp);
    
    Bxx11 = phi_xx(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Bxx12 = x_exp.*x_exp_Dx.*(xn(1,ones(1,P)).^x_exp_Dxx).*(yn(1,ones(1,P)).^y_exp);    
    
    By11 = phi_y(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    By12 = y_exp.*(xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp_Dy);
    
    Byy11 = phi_yy(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Byy12 = y_exp.*y_exp_Dy.*(xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp_Dyy);
    
    Bxy11 = phi_xy(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Bxy12 = x_exp.*y_exp.*(xn(1,ones(1,P)).^x_exp_Dxy).*(yn(1,ones(1,P)).^y_exp_Dxy);
    
    Blap11 = phi_lap(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
    Blap12 = x_exp.*x_exp_Dx.*(xn(1,ones(1,P)).^x_exp_Dxx).*(yn(1,ones(1,P)).^y_exp)+y_exp.*y_exp_Dy.*(xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp_Dyy);
    
    A =   [A11   A12;   A12'   zeros(P)];
    Ax =  [Ax11  Ax12;  Ax12'  zeros(P)];
    Axx = [Axx11 Axx12; Axx12' zeros(P)];
    Ay =  [Ay11  Ay12;  Ay12'  zeros(P)];
    Ayy = [Ayy11 Ayy12; Ayy12' zeros(P)];
    Axy = [Axy11 Axy12; Axy12' zeros(P)];
    
    B = [B11 B12];
    
    format short
    weights = B/A;
    
%%% I believe the issue is in using this method for computing the
%%% derivatives. This method fits an interpolation at each point and then
%%% uses that interpolant to compute the derivative. The issue most likely
%%% arises when computing the derivative of the points directly on the 
%%% y = -x line because there is a cusp/ridge along this line. So, when we
%%% try to resolve derivatives along the ridge using interpolants, there is some
%%% undefined behavior which leads to the solution blowing up along this
%%% ridge.

%%% In my code, I discretize our operators and then use that discretization
%%% at each iteration to compute the new derivatives. 



    del_x  = weights*Ax(1,:)';
    del_y  = weights*Ay(1,:)';
    del_xx = weights*Axx(1,:)';
    del_yy = weights*Ayy(1,:)';
    del_xy = weights*Axy(1,:)';
    
    kappa = abs((del_xx*(del_y^2) - 2*del_x*del_y*del_xy + del_yy*(del_x^2)) / ((del_x^2) + (del_y^2)));
    
    zNew(index) = z(idx_band(index,1)) + dt*kappa;
    
% %     Bx = [Bx11 Bx12];
% %     By = [By11 By12];
% %     Blap = [Blap11 Blap12];
% %     
% %     Dx_local = Bx/A;
% %     Dy_local = By/A;
% %     Dlap_local = Blap/A;
% %     
% %     Dx(index,idx_band(index,:)) = Dx_local(1,1:n);
% %     Dy(index,idx_band(index,:)) = Dy_local(1,1:n);
% %     Dlap(index,idx_band(index,:)) = Dlap_local(1,1:n);
% %     
end
% % % figure(1)
% % % % plots new function values on background nodes
% % % scatter(x,y,[],zNew','filled')
% % % drawnow


z = zNew;



%% CoulNewton iterations to visualize zls
% for iter_coulNewt = 1:4
% 
% [idx_zls,dist_zls] = knnsearch(gamma,gamma,'k',n); %KNN amongst zero level set nodes
% 
% [idx_zlsband,dist_zlsband] = knnsearch([x y],gamma,'k',n); %KNN for band closest to zero level set
% %time
% 
% for index = 1:length(gamma)
% isites = [dsites(idx_zlsband(index,:),1)'; dsites(idx_zlsband(index,:),2)'];
% 
% 
%     xn2 = dsites(idx_zlsband(index,:),1);
%     yn2 = dsites(idx_zlsband(index,:),2);
% 
%     xn2 = xn2 - xn2(1,1).*ones(length(xn2),1);
%     yn2 = yn2 - yn2(1,1).*ones(length(yn2),1);
%     
%     [X1,X2]=meshgrid(xn2);
%     [Y1,Y2]=meshgrid(yn2);
%     
%     A11 = phi(X1,X2,Y1,Y2,m);
%     A12 = (xn2(:,ones(1,P)).^x_exp(ones(1,n),:)).*(yn2(:,ones(1,P)).^y_exp(ones(1,n),:));
%     
%     Ax11 = phi_x(X1,X2,Y1,Y2,m);
%     Ax12 = x_exp.*(xn2(:,ones(1,P)).^x_exp_Dx).*(yn2(:,ones(1,P)).^y_exp);
%     
%     Axx11 = phi_xx(X1,X2,Y1,Y2,m);
%     Axx12 = x_exp.*x_exp_Dx.*(xn2(:,ones(1,P)).^x_exp_Dxx).*(yn2(:,ones(1,P)).^y_exp);
%     
%     Ay11 = phi_y(X1,X2,Y1,Y2,m);
%     Ay12 = y_exp.*(xn2(:,ones(1,P)).^x_exp).*(yn2(:,ones(1,P)).^y_exp_Dy);
%     
%     Ayy11 = phi_yy(X1,X2,Y1,Y2,m);
%     Ayy12 = y_exp.*y_exp_Dy.*(xn2(:,ones(1,P)).^x_exp).*(yn2(:,ones(1,P)).^y_exp_Dyy);
%     
%     Axy11 = phi_xy(X1,X2,Y1,Y2,m);
%     Axy12 = x_exp.*y_exp.*(xn2(:,ones(1,P)).^x_exp_Dxy).*(yn2(:,ones(1,P)).^y_exp_Dxy);
%     
%     
%     B11 = z(idx_zlsband(index,:));
%     B12 = (xn2(1,ones(1,P)).^x_exp).*(yn2(1,ones(1,P)).^y_exp);
%     
%     Bx11 = phi_x(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
%     Bx12 = x_exp.*(xn2(1,ones(1,P)).^x_exp_Dx).*(yn2(1,ones(1,P)).^y_exp);
%     
%     Bxx11 = phi_xx(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
%     Bxx12 = x_exp.*x_exp_Dx.*(xn2(1,ones(1,P)).^x_exp_Dxx).*(yn2(1,ones(1,P)).^y_exp);    
%     
%     By11 = phi_y(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
%     By12 = y_exp.*(xn2(1,ones(1,P)).^x_exp).*(yn2(1,ones(1,P)).^y_exp_Dy);
%     
%     Byy11 = phi_yy(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
%     Byy12 = y_exp.*y_exp_Dy.*(xn2(1,ones(1,P)).^x_exp).*(yn2(1,ones(1,P)).^y_exp_Dyy);
%     
%     Bxy11 = phi_xy(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
%     Bxy12 = x_exp.*y_exp.*(xn2(1,ones(1,P)).^x_exp_Dxy).*(yn2(1,ones(1,P)).^y_exp_Dxy);
%     
%     Blap11 = phi_lap(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
%     Blap12 = x_exp.*x_exp_Dx.*(xn2(1,ones(1,P)).^x_exp_Dxx).*(yn2(1,ones(1,P)).^y_exp)+y_exp.*y_exp_Dy.*(xn2(1,ones(1,P)).^x_exp).*(yn2(1,ones(1,P)).^y_exp_Dyy);
%     
%     A = [A11 A12;A12' zeros(P)];
% %    Ax = [Ax11 Ax12;Ax12' zeros(P)];
% %    Axx = [Axx11 Axx12;Axx12' zeros(P)];
% %    Ay = [Ay11 Ay12;Ay12' zeros(P)];
% %    Ayy = [Ayy11 Ayy12;Ayy12' zeros(P)];
% %    Axy = [Axy11 Axy12;Axy12' zeros(P)];
%     
%     B = [B11 B12];
%     
%     
%     weights = B/A;
%     
%     XYinterp_gamma = gamma(index,:) - [dsites(idx_zlsband(index,1),1) dsites(idx_zlsband(index,1),2)];
%     [Xinterp_gamma1, Xinterp_gamma2] = meshgrid(xn2,XYinterp_gamma(:,1));
%     [Yinterp_gamma1, Yinterp_gamma2] = meshgrid(yn2,XYinterp_gamma(:,2));
%      
%      
%     Ainterp11_gamma = phi(Xinterp_gamma1,Xinterp_gamma2, Yinterp_gamma1,Yinterp_gamma2, m);
%     Ainterp12_gamma = (Xinterp_gamma2(:,ones(1,P)).^x_exp(ones(1,n),:)).*(Yinterp_gamma2(:,ones(1,P)).^y_exp(ones(1,n),:));
%     Ainterp_gamma = [Ainterp11_gamma(1,:)  Ainterp12_gamma(1,:)];
%      
%     Axinterp11_gamma = phi_x(Xinterp_gamma1,Xinterp_gamma2, Yinterp_gamma1,Yinterp_gamma2, m);
%     Axinterp12_gamma = x_exp.*(Xinterp_gamma2(:,ones(1,P)).^x_exp_Dx).*(Yinterp_gamma2(:,ones(1,P)).^y_exp);
%     Axinterp_gamma = [Axinterp11_gamma Axinterp12_gamma];
%     
%     
%     Axxinterp11_gamma = phi_xx(Xinterp_gamma1,Xinterp_gamma2, Yinterp_gamma1,Yinterp_gamma2, m);
%     Axxinterp12_gamma = x_exp.*x_exp_Dx.*(Xinterp_gamma2(:,ones(1,P)).^x_exp_Dxx).*(Yinterp_gamma2(:,ones(1,P)).^y_exp);
%     Axxinterp_gamma = [Axxinterp11_gamma Axxinterp12_gamma];
%     
%     Ayinterp11_gamma = phi_y(Xinterp_gamma1,Xinterp_gamma2, Yinterp_gamma1,Yinterp_gamma2, m);
%     Ayinterp12_gamma = y_exp.*(Xinterp_gamma2(1,ones(1,P)).^x_exp).*(Yinterp_gamma2(1,ones(1,P)).^y_exp_Dy);
%     Ayinterp_gamma = [Ayinterp11_gamma Ayinterp12_gamma];
%     
%     Ayyinterp11_gamma = phi_yy(Xinterp_gamma1,Xinterp_gamma2, Yinterp_gamma1,Yinterp_gamma2, m);
%     Ayyinterp12_gamma = y_exp.*y_exp_Dy.*(Xinterp_gamma2(:,ones(1,P)).^x_exp).*(Yinterp_gamma2(:,ones(1,P)).^y_exp_Dyy);
%     Ayyinterp_gamma = [ Ayyinterp11_gamma Ayyinterp12_gamma ];
%     
%     Axyinterp11_gamma = phi_xy(Xinterp_gamma1,Xinterp_gamma2, Yinterp_gamma1,Yinterp_gamma2, m);
%     Axyinterp12_gamma = x_exp.*y_exp.*(Xinterp_gamma2(:,ones(1,P)).^x_exp_Dxy).*(Yinterp_gamma2(:,ones(1,P)).^y_exp_Dxy);
%     Axyinterp_gamma = [Axyinterp11_gamma Axyinterp12_gamma];
%     
%    
%     
%     
%     func_intp = Ainterp_gamma*weights';
%     x_f_intp = Axinterp_gamma *weights';
%     xx_f_intp = Axxinterp_gamma *weights';
%     y_f_intp = Ayinterp_gamma *weights';
%     yy_f_intp = Ayyinterp_gamma *weights';
%     xy_f_intp = Axyinterp_gamma *weights';
%     
%     
%      
%     eta = 0.8;
% %%%     newton =  (func_intp*((xx_f_intp*(y_f_intp^2)-2*x_f_intp*y_f_intp*xy_f_intp+yy_f_intp*(x_f_intp^2))/((x_f_intp^2+y_f_intp^2)^(2/3))));
%     normgradf = sqrt((x_f_intp).^2+(y_f_intp).^2);
%     scl_xgradf = x_f_intp./normgradf;
%     scl_ygradf = y_f_intp./normgradf;
% 
%     nu = 1;
%     mu = 1;
%     Xcoulomb = 0;
%     Ycoulomb = 0;
%     
%     for nn = 2:n
%          
%         dist = sqrt((gamma(index,1) - gamma(idx_zls(index,nn),1))^2+(gamma(index,2)-gamma(idx_zls(index,nn),2))^2);
%         Xcoulomb = Xcoulomb + (mu*(gamma(index,1) - gamma(idx_zls(index,nn),1))/(dist^nu));
%         Ycoulomb = Ycoulomb + (mu*(gamma(index,2) - gamma(idx_zls(index,nn),2))/(dist^nu));
% 
%     end
% 
%     gammaNew(index,:) = gamma(index,:) + dt.*[Xcoulomb Ycoulomb]-dt.*eta.*func_intp.*[scl_xgradf scl_ygradf];
%    
% 
% end
% 
% 
% gamma = gammaNew;
% end


%% Plot progress of zls over background nodes
framrat = 5;
if mod(time,framrat) == 0
    figure(1)
    hold off

    scatter(x(:),y(:),70,z,'filled')
    
%     hold on
%     plot(gamma(:,1),gamma(:,2),'rx')
    axis equal
    colorbar
    caxis([-0.1,0.3])
    drawnow
    
    time
end

end
