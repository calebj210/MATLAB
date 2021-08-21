% Get curve nodes
close all
clear all

N_Gamma = 300;%Approximate number of nodes
p = 4;% Polynomial Degree;
m = 7;% Order of the PHS --or-- shape parameter of GA
n = (p+2)*(p+1);% Number of neighbors (twice the number of polynomial terms)
P = (p+2)*(p+1)/2;% Number of polynomial terms

%% RBF Definition --- PHS
phi =@(X1,X2,Y1,Y2,alpha) sqrt((X1-X2).^2+(Y1-Y2).^2).^alpha;
phi_x =@(X1,X2,Y1,Y2,alpha) -alpha*(X1-X2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_y =@(X1,X2,Y1,Y2,alpha) -alpha*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-2);
phi_xx =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(X1-X2).^2+(Y1-Y2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
phi_yy =@(X1,X2,Y1,Y2,alpha) alpha*((alpha-1)*(Y1-Y2).^2+(X1-X2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(alpha-4);
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


dR = 0:0.001:0.1;


frames = struct('cdata', cell(1, length(dR)), 'colormap', cell(1, length(dR)));  %preallocate structure with 100 elements
count = 1;


for ind_R=1:length(dR)
    
    x = [R*cos(t)+1;R*cos(t)+3.1-dR(ind_R)];
    y = [R*sin(t)+1;R*sin(t)+1];
    nx = [cos(t);cos(t)];
    ny = [sin(t);sin(t)];
    
    [idx,dist] = knnsearch([x y],[x y],'k',n);
    
    Dx = sparse(N_Gamma);
    Dy = sparse(N_Gamma);
    
    dt = 0.01;
    
    Ne_sqrt = 100;
    Ne = Ne_sqrt^2;
    for index=1%:N_Gamma
        xn = x(idx(index,:));
        yn = y(idx(index,:));
        
        [X1,X2]=meshgrid(xn);
        [Y1,Y2]=meshgrid(yn);
        dx = max(xn)-min(xn);
        dy = max(yn)-min(yn);
        
        xe = linspace(min(xn)-3*dx,max(xn)+3*dx,Ne_sqrt)';
        ye = linspace(min(yn)-3*dy,max(yn)+3*dy,Ne_sqrt)';
        
        [Xe,Ye] = meshgrid(xe,ye);
        xe = Xe(:); ye=Ye(:);
        
        [Xe1,Xe2]=meshgrid(xn,xe);
        [Ye1,Ye2]=meshgrid(yn,ye);
        
        A11 = phi(X1,X2,Y1,Y2,m);
        A12 = (xn(:,ones(1,P)).^x_exp(ones(1,n),:)).*(yn(:,ones(1,P)).^y_exp(ones(1,n),:));
        
        Bx11 = phi_x(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
        Bx12 = x_exp.*(xn(1,ones(1,P)).^x_exp_Dx).*(yn(1,ones(1,P)).^y_exp);
        
        By11 = phi_y(X1(1,:),X2(1,:),Y1(1,:),Y2(1,:),m);
        By12 = y_exp.*(xn(1,ones(1,P)).^x_exp).*(yn(1,ones(1,P)).^y_exp_Dy);
        
        A = [A11 A12;Bx11 Bx12;By11 By12;A12' zeros(P)];
        Bx = [Bx11 Bx12];
        By = [By11 By12];
        
        Ae11 = phi(Xe1,Xe2,Ye1,Ye2,m);
        Ae12 = (xe(:,ones(1,P)).^x_exp(ones(1,Ne),:)).*(ye(:,ones(1,P)).^y_exp(ones(1,Ne),:));
        
        Bxe11 = phi_x(Xe1,Xe2,Ye1,Ye2,m);
        Bxe12 = x_exp.*(xe(:,ones(1,P)).^x_exp_Dx).*(ye(:,ones(1,P)).^y_exp);
        
        Bye11 = phi_y(Xe1,Xe2,Ye1,Ye2,m);
        Bye12 = y_exp.*(xe(:,ones(1,P)).^x_exp).*(ye(:,ones(1,P)).^y_exp_Dy);
        
        
        fn = [ones(n,1);-nx(index);-ny(index);zeros(P,1)];
        lam = A\fn;
        
        fe = [Ae11 Ae12]*lam;
        fxe = [Bxe11 Bxe12]*lam;
        fye = [Bye11 Bye12]*lam;
        
        Dx_local = Bx/A;
        Dy_local = By/A;
        
        
        Dx(index,idx(index,:)) = Dx_local(1,1:n);
        Dy(index,idx(index,:)) = Dy_local(1,1:n);
        
%         hold on
%         plot3(xn,yn,fn(1:n,1),'bo')
%         plot3(xe,ye,fe,'.')
%         plot3(xe,ye,fe+dt*sqrt(fxe.^2+fye.^2),'.')
%         contour(Xe,Ye,reshape(fe,Ne_sqrt,Ne_sqrt),[1 1],'r-')
        contour(Xe,Ye,reshape(fe+dt*sqrt(fxe.^2+fye.^2),Ne_sqrt,Ne_sqrt),[1 1])
        axis equal
        axis([1.9 2.1 0 2])
        frames(ind_R)=getframe(gcf);
%         pause
    end
end

vw = VideoWriter(sprintf('filename%d.avi', 1));  %taking a guess that you intend to modify the filename each time you write a video
      open(vw);
      writeVideo(vw, frames);
      