function display_function_modified(Df,true_f,a,b,nxy,h,bounds,h_index,radius,sing_info)
%%Displays the hypergeometric function aFb(z)
%%Modified version to display real, imaginary, magnitude of Df, and error

lw = 4;

al = a(end)-1;
bet = b(end)-a(end)-1;

x = h*(nxy(1):nxy(2));
y = h*(nxy(3):nxy(4));
[xr,xi] = meshgrid(x,y(end:-1:1)); z = complex(xr,xi);
Z=z(:);

cntr_x=-nxy(1)+1;
cntr_y=nxy(4)+1;

Df = reshape(Df,-nxy(3)+nxy(4)+1,-nxy(1)+nxy(2)+1);

if isempty(sing_info)==0
    z_mod = z; z_mod(abs(z-sing_info(1))<eps)=NaN+1i*NaN;
    true = true_f(z_mod);
else
    true = true_f(z);
end

error=abs(Df-true)./abs(true);
error(error==0)=eps;
error(isinf(error))=NaN;

bx = h*nxy;


error_in = max(abs(error(abs(z)<radius)))
error_out = max(abs(error(abs(z)>radius)))

ax=figure(1)
hold on
subplot(2,2,1) % Plot the real part using mesh
Df_below = Df(cntr_y:end,:);
Df_below(1,:)=conj(Df_below(1,:));
Df_top = Df(1:cntr_y,:);

mesh(xr(cntr_y:end,:),xi(cntr_y:end,:),real(Df_below), 'edgecolor', 'k');
hold on;
mesh(xr(1:cntr_y,:),xi(1:cntr_y,:),real(Df_top), 'edgecolor', 'k');
xlabel('\itx'); ylabel('\ity');
set( gca,'FontSize', 16);
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(1,:));
plot3(xr(cntr_y,:),xi(cntr_y,:),real(Df(cntr_y,:)),'r','LineWidth',lw); % Highlight real axis
plot3(xr(cntr_y,:),xi(cntr_y,:),real(Df_below(1,:)),'r','LineWidth',lw); % Highlight real axis
title(['$$Re(_aF_b)$$'],'interpreter','latex')

subplot(2,2,2) % Plot the imaginary part using mesh
mesh(xr(cntr_y:end,:),xi(cntr_y:end,:),imag(Df_below), 'edgecolor', 'k');
hold on;
mesh(xr(1:cntr_y,:),xi(1:cntr_y,:),imag(Df_top), 'edgecolor', 'k');
xlabel('\itx'); ylabel('\ity');
set( gca,'FontSize', 16);
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(2,:));
plot3(xr(cntr_y,:),xi(cntr_y,:),imag(Df(cntr_y,:)),'r','LineWidth',lw); % Highlight real axis
plot3(xr(cntr_y,:),xi(cntr_y,:),imag(Df_below(1,:)),'r','LineWidth',lw); % Highlight real axis
title(['$$Im(_aF_b)$$'],'interpreter','latex')



subplot(2,2,3) % Plot the magnitude as a surface
% Increase resolution for smoother display

p = surf(xr,xi,abs(Df),angle(-Df)); % Display the surface
set (p,'EdgeColor','none'); colormap hsv(600); hold on;
xlabel('\itx'); ylabel('\ity'); %title('Magnitude, with phase plot')
title(['$$|_aF_b|$$'],'interpreter','latex')
set( gca,'FontSize', 16);
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(3,:));


plot3(xr(cntr_y,:),xi(cntr_y,:),abs(Df(cntr_y,:)),'k','LineWidth',lw); % Highlight real axis
axes('Position',[0.05 0.05 .17 .17]) % Add color wheel for phase information
[th,r] = meshgrid(linspace(-pi,pi),linspace(0,1));
[X,Y] = pol2cart(th+pi,r);

contourf(X,Y,th,100,'linestyle','none'); hold on % Show wheel colors
plot([-1 1],[0,0],'k'); plot([0 0],[-1,1],'k'); % Show Re and Im axes
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k'); colormap hsv(600);
axis equal; axis off

if ~isnan(max(max(error))) %This will happen if (i) there is no true val info or
    %(ii) if we are displaying a different sheet

    ax4=subplot(2,2,4)
    [M,c] = contourf(real(z),imag(z),log10(abs(error)));
    c.LineWidth = 1;
    xlabel('real z','interpreter','latex');
    ylabel('imag z','interpreter','latex');
    set( gca,'FontSize', 16);
    colorbar
    title({['$$Relative \ error$$'];['$$a = $$' num2str(a) ' $$ \ \ b = $$' num2str(b)]},'interpreter','latex')
    hold on
    plot(real(z),imag(z),'k.')
    rectangle('Position',[-radius -radius 2*radius 2*radius],'Curvature',[1 1],'EdgeColor','r','LineWidth',2)
    colormap(ax4,"parula")
    set( gca,'FontSize', 16);
    colormap(ax4,"parula")
    axis equal

else

    subplot(2,2,3)
    colormap hsv(600); hold on;
    xlabel('\itx'); ylabel('\ity'); %title('Magnitude, with phase plot')
    title(['$$|D^{\alpha} f|$$'],'interpreter','latex')
    set( gca,'FontSize', 16);
    xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(3,:));

    axes('Position',[0.05 0.05 .17 .17]) % Add color wheel for phase information
    [th,r] = meshgrid(linspace(-pi,pi),linspace(0,1));
    [X,Y] = pol2cart(th+pi,r);


    contourf(X,Y,th,100,'linestyle','none'); hold on % Show wheel colors
    plot([-1 1],[0,0],'k'); plot([0 0],[-1,1],'k'); % Show Re and Im axes
    plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k'); colormap hsv(600);
    axis equal; axis off


    ax4=subplot(2,2,4)
    colormap(ax4,"parula")
end

ax.Position = [100 100 800 800];
