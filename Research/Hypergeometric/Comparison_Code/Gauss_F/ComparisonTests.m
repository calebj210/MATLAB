%%
% Create method comparison table for 2F1 Hypergeometric functions
%
%
% Author: Caleb Jacobs
% Date last modified: 19-Sep-2022
%%

% Define 30 test cases (a,b,c,z)
test = [0.1,  0.2,  0.3,   0.5; ...
       -0.1,  0.2,  0.3,   0.5; ...
        0.1,  0.2, -0.3,  -0.5 + 0.5i; ...
        1e-8, 1e-8, 1e-8,  1e-6; ...
        1e-8,-1e-6, 1e-12,-1e-10 + 1e-12i; ...
        1,    10,   1,     0.5 + 1e-9i; ...
        1, -1 + 1e-12i, 1, -0.8; ...
        2 + 8i, 3 - 5i, sqrt(2) - pi*1i, 0.75; ...
        100, 200, 350, 1i; ...
        2 + 1e-9, 3, 5, -0.75; ...
        -2, -3, -5 + 1e-9, 0.5; ...
        -1, -1.5, -2 - 1e-15, 0.5; ...
        500, -500, 500, 0.75; ...
        500, 500, 500, 0.75; ...
        -1000, -2000, -4000.1, -0.5; ...
        -100, -200, -300 + 1e-9, 0.5*sqrt(2);
        300, 10, 5, 0.5;
        5, -300, 10, 0.5;
        10, 5, -300.5, 0.5;
        2 + 200i, 5, 10, 0.6;
        2 + 200i, 5 - 100i, 10 + 500i, 0.8;
        2, 5, 10 - 500i, -0.8;
        2.25, 3.75, -0.5, -1;
        1, 2, 4 + 3i, 0.6 - 0.8i;
        1, 0.9, 2, exp(1i * pi/3);
        1, 1, 4, exp(1i * pi/3);
        -1, 0.9, 2, exp(-1i * pi/3);
        4, 1.1, 2, 0.5 + (0.5 * sqrt(3) - 0.01) * 1i;
        5, 2.2, -2.5, 0.49 + 0.5 * sqrt(3) * 1i;
        2/3, 1, 4/3, exp(1i * pi/3)];

% Compute "true" solutions
tru = zeros(30, 1);
for i = 1 : 30
    a = test(i, 1);
    b = test(i, 2);
    c = test(i, 3);
    z = test(i, 4);
    tru(i) = hypergeom([a,b], c, z);
end

% Begin testing
results = zeros(30, 7);
for i = 1 : 30
    a = test(i, 1);
    b = test(i, 2);
    c = test(i, 3);
    z = test(i, 4);

    gam = gamfun(c);
    
    % Pochammer
    val = gam * hypfun_F_pochham(a, b, c, z, 1000);
    results(i, 1) = log10(abs(tru(i) - val));

    % Gauss-Jacobi
    if (real(c - b - 1) >= -1 && real(b - 1) >= -1)
        val = gam * hypfun_F_gjquad(a, b, c, z, 16);
        results(i, 2) = log10(abs(tru(i) - val));
    end
    
    % Taylor (a)
    val = gam * hypfun_F_taylora(a, b, c, z, 1e-15);
    results(i, 3) = log10(abs(tru(i) - val));
    
    % Taylor (b)
    val = gam * hypfun_F_taylorb(a, b, c, z, 1e-15);
    results(i, 4) = log10(abs(tru(i) - val));

    % Single fraction
    val = gam * hypfun_F_singlefraction(a, b, c, z, 1e-15);
    results(i, 5) = log10(abs(tru(i) - val));

    % Buhring
    val = gam * hypfun_F_buhring(a, b, c, z, 0.8, 1e-15);
    results(i, 6) = log10(abs(tru(i) - val));

    % Michel and Stoitsov
    val = gam * hypfun_F_michelstoitsov(a, b, c, z, 1e-15);
    results(i, 7) = log10(abs(tru(i) - val));
end

sprintf("Poc\tG-J\tT-a\tT-b\tSiF\tBuh\tM-S")
results = [[1:30]' ceil(abs(min(results, 0)))];
results = ["Test #", "Poc", "G-J", "T-a", "T-b", "SiF", "Buh", "M-S"; results]