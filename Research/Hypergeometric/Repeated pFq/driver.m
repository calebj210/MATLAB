clear; close all;
clc
format long e

nx = 50;% number of nodes along real axis
h_list = .05;

bounds = [ -2, 2; -2, 2; 0, 3]; % Lower and upper in the three display

for h_index = 1:length(h_list)
    h = h_list(h_index);
    bx = [-nx*h/2,nx*h/2,-nx*h/2,nx*h/2]; % Domain surrounding origin to be displayed
    nxy=sign(bx).*ceil(abs(bx)/h);

    %Approximate pFq([a];[b];x)
    % Define 30 test cases (a,b,c,z)
    test = [0.1,  0.2,  0.3,   0.5; ...                  % 1
       -0.1,  0.2,  0.3,   0.5; ...                      % 2
        0.1,  0.2, -0.3,  -0.5 + 0.5i; ...               % 3
        1e-8, 1e-8, 1e-8,  1e-6; ...                     % 4
        1e-8,-1e-6, 1e-12,-1e-10 + 1e-12i; ...           % 5
        1,    10,   1,     0.5 + 1e-9i; ...              % 6
        1, -1 + 1e-12i, 1, -0.8; ...                     % 7
        2 + 8i, 3 - 5i, sqrt(2) - pi*1i, 0.75; ...       % 8
        100, 200, 350, 1i; ...                           % 9
        2 + 1e-9, 3, 5, -0.75; ...                       % 10
        -2, -3, -5 + 1e-9, 0.5; ...                      % 11
        -1, -1.5, -2 - 1e-15, 0.5; ...                   % 12
        500, -500, 500, 0.75; ...                        % 13
        500, 500, 500, 0.75; ...                         % 14
        -1000, -2000, -4000.1, -0.5; ...                 % 15
        -100, -200, -300 + 1e-9, 0.5*sqrt(2);            % 16
        300, 10, 5, 0.5;                                 % 17
        5, -300, 10, 0.5;                                % 18
        10, 5, -300.5, 0.5;                              % 19
        2 + 200i, 5, 10, 0.6;                            % 20
        2 + 200i, 5 - 100i, 10 + 500i, 0.8;              % 21
        2, 5, 10 - 500i, -0.8;                           % 22
        2.25, 3.75, -0.5, -1;                            % 23
        1, 2, 4 + 3i, 0.6 - 0.8i;                        % 24
        1, 0.9, 2, exp(1i * pi/3);                       % 25
        1, 1, 4, exp(1i * pi/3);                         % 26
        -1, 0.9, 2, exp(-1i * pi/3);                     % 27
        4, 1.1, 2, 0.5 + (0.5 * sqrt(3) - 0.01) * 1i;    % 28
        5, 2.2, -2.5, 0.49 + 0.5 * sqrt(3) * 1i;         % 29
        2/3, 1, 4/3, exp(1i * pi/3)];                    % 30
%     a = [1.6]; %1F1
    
    k = 3;

    a = [test(k,1) test(k,2)]; %2F1
    b = [test(k,3)];

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

    [Df,radius] = compute_pFq_circle(f,a,b,nxy,h,sing_info);

    display_function(Df,true_f,a,b,nxy,h,bounds,h_index,radius,sing_info); % Create the three subplots

end