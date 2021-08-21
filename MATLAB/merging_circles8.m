% circles colliding test case
% Uses 3 layers (separated by a distance of delta) and no derivative 
% constraint in this version (it doesn't mean that we won't want to add
% that in the future).
% It solves S_t = |grad(S)|, so assuming that the right hand side is
% constant in time, S(t+dt) = S(t) + t |grad(S(t))|.
% We know that the (x,y) points for which S(t,x,y)=0 belong to the curve
% that we want to follow at time t. Now we evolve in time and want to find
% the curve once again. So we need the zero level set of 
% S(t) + t |grad(S(t))|. To do so, we use a Newton iteration. 
% Also, for each point, we rotate the coordinate system for it to become
% the origin and its corresponding normal vectoor, the y-axis.

% What can make it better?
% 1) We need to make it more robust - to get away with large dt values
% 2) There are sometimes some glitches at the intersection between the
% circles - it's from the Newton iteration -- the wrong zero level set is
% chosen. That should be fixed with the correct parameter values and using
% spatial nodes refinement.
% 3) Spatial nodes refinement. We want to use the Coul-Newton on segments
% of the zero level set.

close all
clear all
format long e
warning off
N_Gamma = 100;%Approximate number of nodes
p = 2;% Polynomial Degree;
m = 7;% Order of the PHS --or-- shape parameter of GA
n = (p+2)*(p+1);% Number of neighbors (twice the number of polynomial terms)
P = (p+2)*(p+1)/2;% Number of polynomial terms
q = 1;%Number of constrained nodes
tol = 1e-10; %Tolerance for Newtoon iterations
max_iter = 10; %Newton iterations max iterations
eval_dense_grid = 0;
dt = 0.01;
delta = 1*dt;
Final_time = 2;

%% RBF Definition --- PHS
phi =@(X1,X2,Y1,Y2,alpha) sqrt((X1-X2).^2+(Y1-Y2).^2).^alpha;
phi_x =@(X1,X2,Y1,Y2,alpha) -alpha*(X1-X2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_y =@(X1,X2,Y1,Y2,alpha) -alpha*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_xx =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(X1-X2).^2+(Y1-Y2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_xy =@(X1,X2,Y1,Y2,alpha) alpha*(alpha-2)*(X1-X2).*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_yy =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(Y1-Y2).^2+(X1-X2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_xyy =@(X1,X2,Y1,Y2,alpha) alpha*(alpha-2)*(X2-X1).*((X1-X2).^2+(alpha-3)*(Y1-Y2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-6);
phi_yyy =@(X1,X2,Y1,Y2,alpha) alpha*(alpha-2)*(Y2-Y1).*((alpha-1)*(Y1-Y2).^2+3*(X1-X2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-6);
phi_yyyy =@(X1,X2,Y1,Y2,alpha) alpha*(alpha-2)*(3*((Y1-Y2).^2+(X1-X2).^2).^2+6*(alpha-4)*(Y2-Y1).*((Y1-Y2).^2+(X1-X2).^2)+(alpha-4)*(alpha-6)*(Y1-Y2).^4).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-8);
phi_lap =@(X1,X2,Y1,Y2,alpha) alpha^2*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_LS = @(X1,X2,Y1,Y2,alpha) alpha*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-1);

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

t = linspace(0,2*pi,N_Gamma+1)'; t(end)=[];
R = 1;

Cx1 = 1; Cy1 = 1;
Cx2 = 3; Cy2 = 1;

%Gamma
new_x = [R*cos(t)+Cx1;R*cos(t)+Cx2];
new_y = [R*sin(t)+Cy1;R*sin(t)+Cy2];

F_gamma = zeros(length(new_x),1);

%Normals
new_nx = -[cos(t);cos(t)];
new_ny = -[sin(t);sin(t)];

%Off Gamma - layer 1
new_x_1 = [(R+delta)*cos(t)+Cx1;(R+delta)*cos(t)+Cx2];
new_y_1 = [(R+delta)*sin(t)+Cy1;(R+delta)*sin(t)+Cy2];

%Off Gamma - layer 2
new_x_2 = [(R-delta)*cos(t)+Cx1;(R-delta)*cos(t)+Cx2];
new_y_2 = [(R-delta)*sin(t)+Cy1;(R-delta)*sin(t)+Cy2];

F_gamma = [F_gamma;-delta*ones(length(new_x_1),1);delta*ones(length(new_x_1),1)];


Ne_sqrt = 100;
Ne = Ne_sqrt^2;

%True solution
xt = new_x(:); yt = new_y(:);
nxt = -[cos(t);cos(t)]; nyt = -[sin(t);sin(t)];
txt = -nyt; tyt = nxt;

for time = dt:dt:Final_time
    
    %Approx
    x = [new_x(:);new_x_1(:);new_x_2(:)]; y = [new_y(:);new_y_1(:);new_y_2(:)];
    nx = [new_nx(:);new_nx(:);new_nx(:)]; ny = [new_ny(:);new_ny(:);new_ny(:)];
    tx = -ny; ty = nx;
    
    new_nx=[];new_ny=[];
    
    %True solution
    new_xt=[];new_yt=[];
    
    [idx,dist] = knnsearch([x y],[x y],'k',n);
    
    new_x = [];new_y = [];
    new_x_1 = [];new_y_1 = [];
    new_x_2 = [];new_y_2 = [];
    
    
    for index=1:length(x)/3 %3 layers
        if isfinite(x(index)) & isfinite(y(index))
            xn = x(idx(index,:))-x(index);
            yn = y(idx(index,:))-y(index);
            
            Fn = F_gamma(idx(index,:));
            
            % rotates the cluster for the normal to be the y-axis
            Yn = nx(index)*xn+ny(index)*yn;
            Xn = tx(index)*xn+ty(index)*yn;
                        
            
            [X1,X2]=meshgrid(Xn);
            [Y1,Y2]=meshgrid(Yn);
            dx = max(Xn)-min(Xn);
            dy = max(Yn)-min(Yn);
            
            
            A11 = phi(X1,X2,Y1,Y2,m);
            A12 = (Xn(:,ones(1,P)).^x_exp(ones(1,n),:)).*(Yn(:,ones(1,P)).^y_exp(ones(1,n),:));
            
            Bx11 = phi_x(X1(1:q,:),X2(1:q,:),Y1(1:q,:),Y2(1:q,:),m);
            Bx12 = x_exp.*(Xn(1:q,ones(1,P)).^x_exp_Dx).*(Yn(1:q,ones(1,P)).^y_exp);
            Bx = [Bx11 Bx12];
            
            By11 = phi_y(X1(1:q,:),X2(1:q,:),Y1(1:q,:),Y2(1:q,:),m);
            By12 = y_exp.*(Xn(1:q,ones(1,P)).^x_exp).*(Yn(1:q,ones(1,P)).^y_exp_Dy);
            By = [By11 By12];
            
            Byy11 = phi_yy(X1(1:q,:),X2(1:q,:),Y1(1:q,:),Y2(1:q,:),m);
            Byy12 = y_exp.*y_exp_Dy.*(Xn(1:q,ones(1,P)).^x_exp).*(Yn(1:q,ones(1,P)).^y_exp_Dyy);
            Byy = [Byy11 Byy12];
            
            
            %A = [[[A11 A12 ;A12' zeros(P)] By' Byy' Bx']; [[By ; Byy ; Bx] zeros(3*q)]] ;
            A = [A11 A12 ;A12' zeros(P)] ;
            
            %fn = [zeros(n,1);zeros(P,1);ones(q,1);zeros(q,1);zeros(q,1)];
            fn = [Fn;zeros(P,1)];
            lam = A\fn;
            
            S = @(X,Y) [phi(Xn',X,Yn',Y,m) (X(:,ones(1,P)).^x_exp(1,:)).*(Y(:,ones(1,P)).^y_exp(1,:))]*lam;% phi_y(Xn(1:q)',X,Yn(1:q)',Y,m) phi_yy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_x(Xn(1:q)',X,Yn(1:q)',Y,m) ]*lam;
            Sx = @(X,Y) [phi_x(Xn',X,Yn',Y,m) x_exp(1,:).*(X(:,ones(1,P)).^x_exp_Dx(1,:)).*(Y(:,ones(1,P)).^y_exp(1,:))]*lam;% phi_xy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_xyy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_xx(Xn(1:q)',X,Yn(1:q)',Y,m) ]*lam;
            Sy = @(X,Y) [phi_y(Xn',X,Yn',Y,m) y_exp(1,:).*(X(:,ones(1,P)).^x_exp(1,:)).*(Y(:,ones(1,P)).^y_exp_Dy(1,:))]*lam;% phi_yy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_yyy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_xy(Xn(1:q)',X,Yn(1:q)',Y,m) ]*lam;
            %         Sxx = @(X,Y) [phi_xx(xn',X,yn',Y,m) x_exp_Dx(1,:).*x_exp(1,:).*(X(:,ones(1,P)).^x_exp_Dxx(1,:)).*(Y(:,ones(1,P)).^y_exp(1,:))]*lam;
            %         Sxy = @(X,Y) [phi_xy(xn',X,yn',Y,m) x_exp(1,:).*y_exp(1,:).*(X(:,ones(1,P)).^x_exp_Dx(1,:)).*(Y(:,ones(1,P)).^y_exp_Dy(1,:))]*lam;
            Syy = @(X,Y) [phi_yy(xn',X,yn',Y,m) y_exp_Dy(1,:).*y_exp(1,:).*(X(:,ones(1,P)).^x_exp(1,:)).*(Y(:,ones(1,P)).^y_exp_Dyy(1,:))]*lam;%  phi_yyy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_yyyy(Xn(1:q)',X,Yn(1:q)',Y,m) phi_xyy(Xn(1:q)',X,Yn(1:q)',Y,m) ]*lam;
            
            
            
            
            if eval_dense_grid==1 %More info
                
                disp(['S : ' num2str(S(0,0))])
                disp(['Sx : ' num2str(Sx(0,0))])
                disp(['Sy : ' num2str(Sy(0,0))])
                
                
                xe = linspace(Xn(1)-2*dt,Xn(1)+2*dt,Ne_sqrt)';
                ye = linspace(Yn(1)-2*dt,Yn(1)+2*dt,Ne_sqrt)';
                
                [Xe,Ye] = meshgrid(xe,ye);
                xe = Xe(:); ye=Ye(:);
                
                [Xe1,Xe2]=meshgrid(Xn,xe);
                [Ye1,Ye2]=meshgrid(Yn,ye);
                
                Ae11 = phi(Xe1,Xe2,Ye1,Ye2,m);
                Ae12 = (xe(:,ones(1,P)).^x_exp(ones(1,Ne),:)).*(ye(:,ones(1,P)).^y_exp(ones(1,Ne),:));
                
                Bxe11 = phi_x(Xe1,Xe2,Ye1,Ye2,m);
                Bxe12 = x_exp.*(xe(:,ones(1,P)).^x_exp_Dx).*(ye(:,ones(1,P)).^y_exp);
                
                Bye11 = phi_y(Xe1,Xe2,Ye1,Ye2,m);
                Bye12 = y_exp.*(xe(:,ones(1,P)).^x_exp).*(ye(:,ones(1,P)).^y_exp_Dy);
                
                Bxxe11 = phi_xx(Xe1,Xe2,Ye1,Ye2,m);
                Bxye11 = phi_xy(Xe1,Xe2,Ye1,Ye2,m);
                Byye11 = phi_yy(Xe1,Xe2,Ye1,Ye2,m);
                Bxyye11 = phi_xyy(Xe1,Xe2,Ye1,Ye2,m);
                Byyye11 = phi_yyy(Xe1,Xe2,Ye1,Ye2,m);
                
                
                %fep = [Ae11 Ae12 Bye11(:,1:q) Byye11(:,1:q) Bxe11(:,1:q)]*lam;
                %fxep = [Bxe11 Bxe12  Bxye11(:,1:q) Bxyye11(:,1:q) Bxxe11(:,1:q)]*lam;
                %fyep = [Bye11 Bye12  Byye11(:,1:q) Byyye11(:,1:q) Bxye11(:,1:q)]*lam;
                
                fep = [Ae11 Ae12]*lam;
                fxep = [Bxe11 Bxe12]*lam;
                fyep = [Bye11 Bye12]*lam;
                
                fe = S(xe,ye);
                fxe = Sx(xe,ye);
                fye = Sy(xe,ye);
                
                [norm(fep-fe) norm(fxep-fxe) norm(fyep-fye)]
                
                figure(2)
                hold off
                plot3(Xn(:),Yn(:),dt*ones(n,1),'go')
                hold on
                plot3(Xn(1),Yn(1),dt,'ro')
                mesh(Xe,Ye,reshape(fe+dt*fye,Ne_sqrt,Ne_sqrt))
                
                figure(3)
                hold off
                contour(Xe,Ye,reshape(fe+dt*fye,Ne_sqrt,Ne_sqrt),[0,0])
                hold on
                plot(Xn(:),Yn(:),'go')
                plot(Xn(1),Yn(1),'ro')
                plot(Xn(1),Yn(1)-dt,'bo')
                pause
                
            end
            
            mu = 0;
            flag=1;
            iter_num=0;
            x1p = 0;
            y1p = 0;
            while flag==1 %Newton iteration
                
                y1p = mu;
                
                %General case
                mu = mu - (S(x1p,y1p)+dt*abs(Sy(x1p,y1p)))/(Sy(x1p,y1p)+dt*sign(Sy(x1p,y1p))*Syy(x1p,y1p));
                %True if Sy positive
                %mu = mu - (S(x1p,y1p)+dt*Sy(x1p,y1p))/(Sy(x1p,y1p)+dt*Syy(x1p,y1p));
                
                y1p = mu;
                
                %General case
                fval = S(x1p,y1p)+dt*abs(Sy(x1p,y1p));
                %True if Sy positive
                %fval = S(x1p,y1p)+dt*Sy(x1p,y1p);
                
                
                iter_num = iter_num+1;
                if abs(fval-0)<tol
                    
                    
                    x1p = y1p*ty(index)/(nx(index)*ty(index)-ny(index)*tx(index))+x(index); %Nodes on the curve
                    y1p = -y1p*tx(index)/(nx(index)*ty(index)-ny(index)*tx(index))+y(index);
                    
                    x1p_1 = x1p-delta*nx(index); %nodes off the curve
                    y1p_1 = y1p-delta*ny(index);
                    
                    x1p_2 = x1p+delta*nx(index); %nodes off the curve
                    y1p_2 = y1p+delta*ny(index);
                    
                    
                    new_x = [new_x x1p]; %Nodes on the curve
                    new_y = [new_y y1p];
                    
                    new_x_1 = [new_x_1 x1p_1]; %nodes off the curve
                    new_y_1 = [new_y_1 y1p_1];
                    
                    new_x_2 = [new_x_2 x1p_2]; %nodes off the curve
                    new_y_2 = [new_y_2 y1p_2];
                    
                    
                    new_nx = [new_nx nx(index)];
                    new_ny = [new_ny ny(index)];
                    
                    flag = 0;
                    
                elseif iter_num == max_iter
                    
                    flag = 0;
                    
                    new_x = [new_x NaN];
                    new_y = [new_y NaN];
                    
                    new_x_1 = [new_x_1 NaN];
                    new_y_1 = [new_y_1 NaN];
                    
                    new_x_2 = [new_x_2 NaN];
                    new_y_2 = [new_y_2 NaN];
                    
                    new_nx = [new_nx NaN];
                    new_ny = [new_ny NaN];
                    
                end
            end
        else
            new_x = [new_x NaN];
            new_y = [new_y NaN];
            
            new_x_1 = [new_x_1 NaN];
            new_y_1 = [new_y_1 NaN];
            
            new_x_2 = [new_x_2 NaN];
            new_y_2 = [new_y_2 NaN];
            
            new_nx = [new_nx NaN];
            new_ny = [new_ny NaN];
        end
    end
    
    %True solution
    for index=1:length(xt)
        x1pt = (-time)*tyt(index)/(nxt(index)*tyt(index)-nyt(index)*txt(index))+xt(index); %Going back to our original coord system
        y1pt = -(-time)*txt(index)/(nxt(index)*tyt(index)-nyt(index)*txt(index))+yt(index);
        
        new_xt = [new_xt x1pt];
        new_yt = [new_yt y1pt];
    end
    
    
    error_linf = max(abs(sqrt((new_x-new_xt).^2+(new_y-new_yt).^2)));
    
    disp(['time = ' num2str(time) '  Error (l_inf) = ' num2str(error_linf) ' number of nodes = ' num2str(sum(isfinite(x)))])
    
    figure(4)
    %hold off
    plot(new_x,new_y,'b.')
    axis equal
end



figure(5)
hold off
plot(new_x,new_y,'b.')
hold on
plot(new_xt,new_yt,'o')
