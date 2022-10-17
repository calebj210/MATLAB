function approx = hyp_integrate(za,zb,f,al,be)
%Linked to my CP221006 handout.
%Function that evaluates using end correction the integral
%
% I = int_za^zb (t-za)^al (zb-t)^be f(t) dt
%
%1)f(t) can have poles but otherwise cannot have branch points (this could
% be fixed though)
%2)The integration is performed along a straight line from za to zb
%3)The end correction stencil is a circle centered at za (zb) with a radius
% ra (rb) and n nodes around them.
%4)The trapezoidal rule starts from the pth node (p = 1,2,...) and stops at
%the pth node before zb. If p=1, we recover the original trapezoidal rule.
%Since p can be different on either end, we have pa by za and pb by zb.
%5)Parameters to study:
%     1) pa and pb
%     2) ra and rb
%     3) n (the number of nodes around the correction circle)

plot_the_stencil = 0;

%% Parameters
N=15; %number of points inside the domain of integration

h = (zb-za)/(N+1);
ra = 3;%round(2+1/abs(5*h));
rb = 3;%round(2+1/abs(5*h));
n = 20;
%% End of Parameters

fa = @(Z) f(Z).*(Z-zb).^be;
fb = @(Z) f(Z).*(Z-za).^al;
fab = @(Z) f(Z).*((Z-za).^al).*((Z-zb).^be);

t=linspace(0,2*pi,n+1)'; t(end)=[];
zsa = za+abs(h)*ra*exp(1i*t);
zsb = zb+abs(h)*rb*exp(1i*t);

pa=3;
pb=3;

%Correction stencils
Wa = stencil_circ(n,al,h/abs(h),pa,(zsa-za)/(abs(h)));
Wb = stencil_circ(n,be,-h/abs(h),pb,(zsb-zb)/(abs(h)));

zk = za+h*((pa):(N-(pb-1))); %TR path nodes

TR=h*sum(fab(zk));
fa_data = fa(zsa(:));
fb_data = fb(zsb(:));

if plot_the_stencil==1
    plot(real(za),imag(za),'rx'); hold on
    plot(real(zb),imag(zb),'rx')
    plot(real(zk(:)),imag(zk(:)),'.')
    plot(zsa(:),'o')
    plot(zsb(:),'o')
    axis equal
end

%%Correction at the end points. Moving the branch cut when the values
%%function values associated to the stencil are affected by the branch
%%cuts of (z-zb)^be and (z-za)^al

%if za and zb are within each other's circle, then use Taylor expansion
if abs(imag(za-zb))<=abs(h)*max(ra,rb)
    if real(za)<=real(zb)
        if imag(zb)<=imag(za)
            fa_data(imag(zsa(:))<imag(zb))=fa_data(imag(zsa(:))<imag(zb))*exp(2*pi*1i*be);
        else
            fa_data(imag(zsa(:))>=imag(zb))=fa_data(imag(zsa(:))>=imag(zb))*exp(-2*pi*1i*be);
        end
    else
        if imag(za)<=imag(zb)
            fb_data(imag(zsb(:))<imag(za))=fb_data(imag(zsb(:))<imag(za))*exp(2*pi*1i*al);
        else
            fb_data(imag(zsb(:))>=imag(za))=fb_data(imag(zsb(:))>=imag(za))*exp(-2*pi*1i*al);
        end
    end
end

if angle(zb)>0
    Cb = exp(1i*pi*be);
else
    Cb = exp(-1i*pi*be);
end

corr_a = -sum(Wa(:).*(fa_data))*h^(al+1);
corr_b = sum(Wb(:).*(fb_data))*(-h)^(be+1);
approx = Cb*(TR + corr_a + corr_b);
