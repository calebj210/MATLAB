function display_function(f,bx,nx,ny,bounds)
% The input parameters described in the script that calls this routine
lw = 4; % Set LineWidth for highlighting values along real axis
x = linspace(bx(1),bx(2),nx); y = linspace(bx(3),bx(4),ny);
[xr,xi] = meshgrid(x,y); z = complex(xr,xi);

figure(1) % Plot the real part using mesh
mesh(xr,xi,real(f(z))); colormap([0 0 0]); hold on;
xlabel('\itx'); ylabel('\ity'); title('Real part')
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(1,:));
plot3(x,zeros(size(x)),real(f(x)),'r','LineWidth',lw); % Highlight real axis

figure(2) % Plot the imaginary part using mesh
mesh(xr,xi,imag(f(z))); colormap([0 0 0]); hold on;
xlabel('\itx'); ylabel('\ity'); title('Imaginary part')
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(2,:));
plot3(x,zeros(size(x)),imag(f(x)),'r','LineWidth',lw); % Highlight real axis

figure(3) % Plot the magnitude as a surface
% Increase resolution for smoother display
x = linspace(bx(1),bx(2),nx*4); y = linspace(bx(3),bx(4),ny*4);
[xr,xi] = meshgrid(x,y); z = complex(xr,xi);
p = surf(xr,xi,abs(f(z)),angle(-f(z))); % Display the surface
set (p,'EdgeColor','none'); colormap hsv(600); hold on;
xlabel('\itx'); ylabel('\ity'); title('Magnitude, with phase plot')
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(3,:));
plot3(x,zeros(size(x)),abs(f(x)),'k','LineWidth',lw); % Highlight real axis
axes('Position',[0.05 0.05 .17 .17]) % Add color wheel for phase information
[th,r] = meshgrid(linspace(-pi,pi),linspace(0,1));
[X,Y] = pol2cart(th+pi,r);
contourf(X,Y,th,100,'linestyle','none'); hold on % Show wheel colors
plot([-1 1],[0,0],'k'); plot([0 0],[-1,1],'k'); % Show Re and Im axes
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k'); colormap hsv(600);
axis equal; axis off

% Test code for displaying an analytic function; here exp(z)
clear; close all;
bx = [-4,4,-10,10]; % Domain (box) surrounding origin to be displayed
nx = 41; % Number of nodes along real axis for re & im plots
ny = 41; % Number of nodes along imag axis for re & im plots
bounds = [ -60,60; -60,60; 0,60]; % Lower and upper in the three displays
f = @(z) z.^2; % Function to display
display_function(f,bx,nx,ny,bounds); % Create the three subplots


f[z_] = Sin[z];

(*x and y bounds*)
{x1, x2, y1, y2} = {-2, 2, -2, 2}; 

(*Real Curves*)
{re, im, abs} = 
  ReleaseHold[Hold[ParametricPlot3D[{x, 0, a[f[x]]}, {x, x1, x2}, 
      PlotStyle -> Directive[Thickness[0.01], b]]] /. 
    {{a -> Re, b -> Red}, 
     {a -> Im, b -> Red}, 
     {a -> Abs, b -> Black}}]; 

(*Real*)
Show[Plot3D[Re[f[x + I*y]], {x, x1, x2}, {y, y1, y2},  
  AxesLabel -> {"Re(z)", "Im(z)"}, PlotLabel -> Re[f[z]], 
  ImageSize -> Medium, ClippingStyle -> Opacity[.75]], re]

(*Imaginary*)
Show[Plot3D[Im[f[x + I*y]], {x, x1, x2}, {y, y1, y2}, 
  AxesLabel -> {"Re(z)", "Im(z)"}, PlotLabel -> Im[f[z]], 
  ImageSize -> Medium, ClippingStyle -> Opacity[.75]], im]

(*Magnitude and Argument*)
Show[
 Plot3D[Abs[f[x + I*y]], {x, x1, x2}, {y, y1, y2},
  ColorFunction -> Function[{x, y, z}, Hue[Arg[f[x + I*y]]/(2*Pi)]],
  ColorFunctionScaling -> False,
  ClippingStyle -> Opacity[.75], 
  PlotPoints -> 100, ImageSize -> Medium, 
     Epilog -> Inset[ParametricPlot[{u*Cos[t], u*Sin[t]}, 
        {u, 0, 1}, {t, 0, 2*Pi}, 
       ColorFunction -> Function[{x, y, u, t}, Hue[t/(2*Pi)]], 
        ColorFunctionScaling -> False, 
        Axes -> False, BoundaryStyle -> None, 
           Mesh -> None, Frame -> None, ImageSize -> 50], {0.05, 
     0.1}], 
        AxesLabel -> {"Re(z)", "Im(z)"}, PlotLabel -> Abs[f[z]]], abs]