clear; close all;
clc
format long e

nx = 50;% number of nodes along real axis
h_list = .05;

bounds = [ -5,5; -5,5; 0,5]; % Lower and upper in the three display

for h_index = 1:length(h_list)
    h = h_list(h_index);
    bx = [-nx*h/2,nx*h/2,-nx*h/2,nx*h/2]; % Domain surrounding origin to be displayed
    nxy=sign(bx).*ceil(abs(bx)/h);

    %Approximate pFq([a];[b];x)
    a = [1.6]; %1F1
    %a = [2.5 3.6]; %2F1
    b = [1.8];

    if length(a)==1 & length(b)==1 %1F1(a;b;x)
        f = @(x) exp(x);
        sing_info = []; %no singularity
    elseif length(a)==2 & length(b)==1 % 2F1(a;b;x)
        f = @(x) (1-x).^(-a(1));
        sing_info = [1, -a(1)]; %branch point at z=1
    end
    true_f = @(x) hypergeom(a, b, x); %True value
    al = a(end)-1;
    bet = b(end)-a(end)-1;

    [Df,radius] = compute_pFq_circle(f,a,b,nxy,h,sing_info)

%     display_function(Df,true_f,a,b,nxy,h,bounds,h_index,radius,sing_info); % Create the three subplots

end