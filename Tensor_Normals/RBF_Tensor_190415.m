%% Calculating Approximate Normals
%The following code calculates the approximate normals to be used with
%RBF's to solve PDE's on arbitrary surfaces.

%% General Start-up
close all
clear all

%% Control variabls
%if fig = 1  DISPLAY GLOBAL FIGURES (normals at all points)
%if fig = 2 DISPLAY LOCAL FIGURES (nearest neighbors and normals)
%if fig = 3 Display error plots
%if fig = 0  NO FIGURES
% fig = 1;
fig = 1;

% if kond = 1 COMPUTE CONDITIONING
% if kond = 0 NO CONDITIONING COMPUTATION
kond = 0;

%% Selecting Node Distribution
%equispaced node distribution
%nodedist = [200; 320; 600; 660; 920; 1320; 1440; 1940; 2080];%Vector with maxgN numbers of node distributions
%nodedist = [660; 920; 1320; 1440; 1940; 2080; 10080];%Vector with maxgN numbers of node distributions
%nodedist = [10080];%Vector with maxgN numbers of node distributions

%nodedist  = [100; 200; 500; 1000; 2000; 5000; 10000];
%nodedist = [100];

nodedist  = [10000];


%% Select minimum and maximum degree of polynomial to append
% degree of polynomial is directly linked to the number of nearest
% neighbors
% Number of neighbors = twice number of terms in appended polynomial
minpolydeg = 2;
maxpolydeg = 2;
% 
%% Curve Parameterization -Unit Circle
% %---Parameterized Unit Circle Curve ---
% xgamma = @(t) cos(t);
% ygamma = @(t) sin(t);
% xpgamma =@(t) -sin(t);
% ypgamma =@(t) cos(t);
% 
% %---True normal and tangent directions
%   xtrunorms=@(t) cos(t);
%   ytrunorms=@(t) sin(t);
%   xtrutangs = @(t) -sin(t);
%   ytrutangs = @(t) cos(t);
%  
% Curve Parameterization - Pi curve
xgamma=@(t) ((17/31).*sin((235/57)+(-32).*t)+(19/17).*sin((192/55)+(-30).*t)+(47/32).*sin((69/25)+(-29).*t)+(35/26).*sin((75/34)+(-27).*t)+(6/31).*sin((23/10)+(-26).*t)+(35/43).*sin((10/33)+(-25).*t)+(126/43).*sin((421/158)+(-24).*t)+(143/57).*sin((35/22)+(-22).*t)+(106/27).*sin((84/29)+(-21).*t)+(88/25).*sin((23/27)+(-20).*t)+(74/27).*sin((53/22)+(-19).*t)+(44/53).*sin((117/25)+(-18).*t)+(126/25).*sin((88/49)+(-17).*t)+(79/11).*sin((43/26)+(-16).*t)+(43/12).*sin((41/17)+(-15).*t)+(47/27).*sin((244/81)+(-14).*t)+(8/5).*sin((79/19)+(-13).*t)+(373/46).*sin((109/38)+(-12).*t)+(1200/31).*sin((133/74)+(-11).*t)+(67/24).*sin((157/61)+(-10).*t)+(583/28).*sin((13/8)+(-8).*t)+(772/35).*sin((59/16)+(-7).*t)+(3705/46).*sin((117/50)+(-6).*t)+(862/13).*sin((19/8)+(-5).*t)+(6555/34).*sin((157/78)+(-3).*t)+(6949/13).*sin((83/27)+(-1).*t)+(-6805/54).*sin((1/145)+2.*t)+(-5207/37).*sin((49/74)+4.*t)+(-1811/58).*sin((55/43)+9.*t)+(-63/20).*sin((2/23)+23.*t)+(-266/177).*sin((13/18)+28.*t)+(-2/21).*sin((7/16)+31.*t))/1000;
ygamma=@(t) ((70/37).*sin((65/32)+(-32).*t)+(11/12).*sin((98/41)+(-31).*t)+(26/29).*sin((35/12)+(-30).*t)+(54/41).*sin((18/7)+(-29).*t)+(177/71).*sin((51/19)+(-27).*t)+(52/33).*sin((133/52)+(-26).*t)+(59/34).*sin((125/33)+(-26).*t)+(98/29).*sin((18/11)+(-25).*t)+(302/75).*sin((59/22)+(-24).*t)+(104/9).*sin((118/45)+(-22).*t)+(52/33).*sin((133/52)+(-21).*t)+(37/45).*sin((61/14)+(-20).*t)+(143/46).*sin((144/41)+(-19).*t)+(254/47).*sin((19/52)+(-18).*t)+(246/35).*sin((92/25)+(-17).*t)+(722/111).*sin((176/67)+(-16).*t)+(136/23).*sin((3/19)+(-15).*t)+(273/25).*sin((32/21)+(-13).*t)+(229/33).*sin((117/28)+(-12).*t)+(19/4).*sin((43/11)+(-11).*t)+(135/8).*sin((23/10)+(-10).*t)+(205/6).*sin((33/23)+(-8).*t)+(679/45).*sin((55/12)+(-7).*t)+(101/8).*sin((11/12)+(-6).*t)+(2760/59).*sin((40/11)+(-5).*t)+(1207/18).*sin((21/23)+(-4).*t)+(8566/27).*sin((39/28)+(-3).*t)+(12334/29).*sin((47/37)+(-2).*t)+(15410/39).*sin((185/41)+(-1).*t)+(-596/17).*sin((3/26)+9.*t)+(-247/28).*sin((25/21)+14.*t)+(-458/131).*sin((21/37)+23.*t)+(-41/36).*sin((7/8)+28.*t))/1000;
xpgamma =@(t) (-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t);
ypgamma =@(t) (-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t);


%---True normal and tangent directions
  xtrunorms=@(t) (1/2000).*((-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t)).*((1/1000000).*((-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t)).^2+(1/1000000).*((-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t)).^2).^(-1);
  ytrunorms=@(t) (-1/2000).*((1/1000000).*((-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t)).^2+(1/1000000).*((-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t)).^2).^(-1).*((-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t));
  xtrutangs = @(t) (1/1000).*((-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t));
  ytrutangs = @(t)(1/1000).*((-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t));

%% Error Data Vectors
%columns: data for line associated with same number of neighbors
%rows: [1:rows/2]  total number of nodes  used, row number corresponds to
%nodedist number entry mesh
%rows: [rows/2 + 1:end] associated error value
errtan = zeros(2*size(nodedist,1), maxpolydeg-minpolydeg+1); 
errnormx = zeros(2*size(nodedist,1),maxpolydeg-minpolydeg+1);
errnormy = zeros(2*size(nodedist,1),maxpolydeg-minpolydeg+1);
errnorml2 = zeros(2*size(nodedist,1),maxpolydeg-minpolydeg+1);

%% Loop over Node distributions
for nodes = 1:size(nodedist,1)  %loop over node distributions 
    %loading first node distribution
    %imported file gives the following
    %"gamma" = cartesian coordinates of points on unit circle
    %"x" = x-coordinates of background nodes
    %"y" = y-coordinates of background nodes

    %% Import Node Distribution
    %Import distributions from equisapced_nodes
    %Imports "gamma", "X", "Y"
    fileName = ['nodes_cN' num2str(nodedist(nodes))];
%    load (fullfile('Equispaced_nodes_circle',[fileName '.mat'])); %unit circle
    load (fullfile('Equispaced_nodes_pi',[fileName '.mat']));   %pi curve

    %% Initialize Vectors - Take one less of gamma to avoid node overlap
    trunorms = zeros(length(gamma)-1,2);
    trutangs = zeros(length(gamma)-1,2);
    
    %% Loop over number of neighbors by increasing degree of appended polynomial
    for polydeg = minpolydeg:maxpolydeg  %loop over adding neighbors by adding polynomial terms
        %% RBF Initialization constants
        N = length(gamma(:,1))-1;            % Get number of nodes on unit circle from gamma\                   % Degree of appended polynomial
        p = polydeg;                      % Degree of appended polynomial
        m = 3;                             % GA Shape parameter or Order of the PHS 
        %n = (p+2)*(p+1);   % Number of neighbors for 2D RBF (twice the number of polynomial terms)
        n = 2*p+2;          %Number of neightbors for 1D RBF
        %P = (p+2)*(p+1)/2; % Number of polynomial terms in 2D RBF
        P= p + 1;                %Number of polynomial terms in 1D RBF
        %% Initiate storage vectors
        calcnorms = zeros(N,2);%calculated normal direction storage = [x-component, y-component]
        trunorms = zeros(N,2);
        calctangs = zeros(N,2);
        trutangs = zeros(N,2);
        
        %% Condition number storage 
        if kond == 1
        CA = zeros(N,1); %Matrix for storing condition number of A
        CA11 = zeros(N,1); %Compute conditioning for RBF component of A
        CA12 = zeros(N,1); %Compute conditioning for Appended polynomial component of A
        end
        
        %% RBF Definitions ---------------------------------------------------------
        %%%Gaussian%%%
        % phi = @(X1,X2,Y1,Y2,m) exp(-m^2*(((X1-X2).^2)+((Y1-Y2).^2)));
        % phi_x = @(X1,X2,Y1,Y2,m) exp(-m^2*(((X1-X2).^2)+((Y1-Y2).^2)))*(-2*m)...
        %     .*(X1-X2);
        % phi_y = @(X1,X2,Y1,Y2,m) exp(-m^2*(((X1-X2).^2)+((Y1-Y2).^2)))*(-2*m)...
        %     .*(Y1-Y2);
        %%%PHS 2D%%%
        phi =@(X1,X2,Y1,Y2,m) sqrt((X1-X2).^2+(Y1-Y2).^2).^m;
        phi_x =@(X1,X2,Y1,Y2,m) -m*(X1-X2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(m-2);
        phi_y =@(X1,X2,Y1,Y2,m) -m*(Y1-Y2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(m-2);
        phi_xx =@(X1,X2,Y1,Y2,m) m*((m-1)*(X1-X2).^2+(Y1-Y2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(m-4);
        phi_yy =@(X1,X2,Y1,Y2,m) m*((m-1)*(Y1-Y2).^2+(X1-X2).^2).*sqrt((X1-X2).^2+(Y1-Y2).^2).^(m-4);
        phi_lap =@(X1,X2,Y1,Y2,m) m^2*sqrt((X1-X2).^2+(Y1-Y2).^2).^(m-2);
        
        %%%PHS 1D%%%
%         phi1D =@(X1,X2,m) sqrt((X1-X2).^2).^m;
       phi1D =@(X1,X2,m) (abs(X1-X2)).^m;
       phi1D_Dx = @(X1,X2,m) -m.*((abs(X1-X2)).^(m-1));

        


        %% Polynomial exponents 2D
        %Sets up help matrices - polynomial augmenting of radial functions
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
        
        %% Polynomial Exponents 1D
        x1D_exp = [];
        for i_exp = 0:p
            x1D_exp=[x1D_exp i_exp];
        end
        x1D_exp_Dx = x1D_exp-1;
        x1D_exp_Dx(x1D_exp_Dx<0)=0;
        x1D_exp_Dxx = x1D_exp-2;
        x1D_exp_Dxx(x1D_exp_Dxx<0)=0;
        
        %% Constructing Approximate Normals
        X = zeros(3*N,3);
        X(1:N,1:2) = gamma(1:N,:);

        initnorms = zeros(N,2);  %storage for initial normal guess, for background nodes

%computing initnorms
%1)choose point on gamma, and two nearest neighbors on gamma.
%2)connect nearest with a line segment and compute midpoint
%3)construct perpendicular vector to this line segment
%4)normalize this vector and scale to have lenth of line segment connecting
%nearest neighbors.
        for i = 1:(N)
           if i >= 2 && i <= N-1
       
               midpoint = [ ((X(i+1,1)+X(i-1,1))/2) ((X(i+1,2)+X(i-1,2))/2)]; %compute midpoint
               %slope = rise/run
               rise = X(i+1,2)-X(i-1,2);
               run = X(i+1,1)-X(i-1,1);
               %perpendicular slope = -run/rise
               initnorms(i,:) = [ (midpoint(1,1)+ (-1)*rise) (midpoint(1,2)+run)];
               initnorms(i,:) = initnorms(i,:) - midpoint;
               scl = 1/sqrt(initnorms(i,1)^2+initnorms(i,2)^2);  %normalize scalar
               initnorms(i,:) = scl.*initnorms(i,:); %scale vector
               ep = sqrt((midpoint(1,1) - X(i-1,1))^2 + (midpoint(1,2) - X(i-1,2))^2); %compute half length of line segment
               initnorms(i,:) = ep.*initnorms(i,:); %scale init norms
       
           elseif i == 1
       
               midpoint = [ ((X(i+1,1)+X(N,1))/2) ((X(i+1,2)+X(N,2))/2)]; %compute midpoint
               %slope = rise/run
               rise = X(i+1,2)-X(N,2);
               run = X(i+1,1)-X(N,1);
               %perpendicular slope = -run/rise
               initnorms(i,:) = [ (midpoint(1,1)+ (-1)*rise) (midpoint(1,2)+run)];
               initnorms(i,:) = initnorms(i,:) - midpoint;
               scl = 1/sqrt(initnorms(i,1)^2+initnorms(i,2)^2);  %normalize scalar
               initnorms(i,:) = scl.*initnorms(i,:); %scale vector
               ep = sqrt((midpoint(1,1) - X(N,1))^2 + (midpoint(1,2) - X(N,2))^2); %compute half length of line segment
               initnorms(i,:) = ep.*initnorms(i,:); %scale init norms
               
           elseif i == N
       
               midpoint = [ ((X(1,1)+X(i-1,1))/2) ((X(1,2)+X(i-1,2))/2)];%compute midpoint
                %slope = rise/run
               rise = X(1,2)-X(i-1,2);
               run = X(1,1)-X(i-1,1);
               %perpendicular slope = -run/rise
               initnorms(i,:) = [ (midpoint(1,1)+ (-1)*rise) (midpoint(1,2)+run)];
               initnorms(i,:) = initnorms(i,:) - midpoint;
               scl = 1/sqrt(initnorms(i,1)^2+initnorms(i,2)^2);%normalize scalar
               initnorms(i,:) = scl.*initnorms(i,:);  %scale vector
               ep = sqrt((midpoint(1,1) - X(i-1,1))^2 + (midpoint(1,2) - X(i-1,2))^2); %compute half length of line segment
               initnorms(i,:) = ep.*initnorms(i,:);  %scale init norms

   end
end

        %Get initnorms to all point in same direction
        %1) At first point on gamma select initnorm vector 1
        %2) At next point on gamma select initnorm vector 2
        %3) Compute sign((initnorm1).(initnorm2))
        %4) Scale initnorm2 by this value
        %5) Repeat over all nodes
		for i = 1:2N;
            if  i<N || N+1<= i <= 2N-1;
              dir = sign((initnorms(i,1)*initnorms(i+1,1) + initnorms(i,2)*initnorms(i+1,2)));
                initnorms(i+1,:) = dir.*initnorms(i+1,:);        
            elseif i == N || i == 2N;
                dir = sign((initnorms(i,1)*initnorms(1,1) + initnorms(i,2)*initnorms(1,2)));
                initnorms(1,1) = dir.*initnorms(1,1);  
                initnorms(1,2) = dir.*initnorms(1,2);
            end
        end

%         %first vector points inwards
%         for i = 1:N
%             %Background nodes inside gamma
%             X(N+i,1:2) = X(i,1:2) + initnorms(i,:);
%             X(N+i,3) = -1;
%             %Background nodes outside gamma
%             X(2*N+i,1:2) = X(i,1:2) - initnorms(i,:);
%             X(2*N+i,3) = +1;
%         end

        %Plot of function values at gamma nodes and background nodes
%         if fig == 1
%         figure
%          title('X nodes')
%          axis equal
%          hold on
%          scatter(X(:,1),X(:,2),25,X(:,3),'filled')
%          
%          figure
%          title('init norms')
%          axis equal
%          hold on
%          delta = 5;
%          plot((X(1:N,1)),(X(1:N,2)),'ko')
%          for i = 1:N
%          plot([X(i,1) (X(i,1)-(delta*initnorms(i,1)))],[X(i,2) (X(i,2)-(delta*initnorms(i,2)))],'r')
%          end
%         end


        % % %uncomment to plot initnorms point by point
        % figure
        % axis equal
        % for i = 1:N
        % hold on
        % plot((X(1:N,1)),(X(1:N,2)),'w*')
        % plot((X(1:N,1)),(X(1:N,2)),'ko')
        % plot((X(i,1)+initnorms(i,1)),(X(i,2)+initnorms(i,2)),'r*')
        % plot((X(i,1)),(X2(i,2)),'g*')
        % if i >1
        % for j = 1:i-1
        %     plot((X(j,1)),(X(j,2)),'mx')
        % end
        % end
        % pause
        % end
        
%normalize normals
           scal = 1./sqrt(initnorms(:,1).^2 + initnorms(:,2).^2);
           initnorms(:,1) = -scal.*initnorms(:,1);
           initnorms(:,2) = -scal.*initnorms(:,2);
           
 

        
        %% Compute initial tangent vectors
        rotang = -pi/2;
        inittangs = zeros(N,2);
        %Apply rotation matrix to each component of initnorms
        for i = 1:N
            inittangs(i,1) = initnorms(i,1)*cos(rotang)-initnorms(i,2)*sin(rotang);
            inittangs(i,2) = initnorms(i,1)*sin(rotang)+initnorms(i,2)*cos(rotang);
        end
        
%         if fig ==1
%             figure
%             title("init tangs")
%             hold on
%             axis equal
%             delta = 0.5;
%             plot((X(1:N,1)),(X(1:N,2)),'ko')
%             for i = 1:N
%             plot([X(i,1) (X(i,1)+(delta*inittangs(i,1)))],[X(i,2) (X(i,2)+(delta*inittangs(i,2)))],'b')
%             end
%         end
            
        %% Compute Nearest Neighbors
%         dsites = [X(:,1) X(:,2)];
%         Nn = [X(:,1) X(:,2)];
        dsites = [X(1:N,1) X(1:N,2)];
        Nn = [X(1:N,1) X(1:N,2)];
       [idx,dist] = knnsearch(X(1:N,[1 2]),gamma,'k',n);

% 
%uncomment to plot nearest neighbors
% if fig ==1
% figure
% title('nearest neighbors check')
% axis equal
% for i = 1:N
% hold on
% plot((X(1:N,1)),(X(1:N,2)),'w*')
% plot((X(1:N,1)),(X(1:N,2)),'ko')
% plot((X(N+1:length(X(:,1)),1)),(X(N+1:length(X(:,1)),2)),'b*')
% plot((X(idx(i,:),1)),(X(idx(i,:),2)),'r*')
% plot((X(i,1)),(X(i,2)),'g*')
% if i >1
% for j = 1:i-1
%     plot((X(j,1)),(X(j,2)),'mx')
% end
% end
% pause
% end
% end

%% Populating RBF Interpolation Matrix
for index = 1:N
    %select nearest neightbor nodes       
    isites = [ (dsites(idx(index,:),1))' ; (dsites(idx(index,:),2))'];
    %recenter nodes at the origin    
    isites(1,:) = isites(1,:) - (isites(1,1)*ones(1,length(isites)));
    isites(2,:) = isites(2,:) - (isites(2,1)*ones(1,length(isites)));
    %Form change of basis matrix
    vardet = inittangs(index,1)*initnorms(index,2)-inittangs(index,2)*initnorms(index,1);
    changemat = (1/vardet).*[initnorms(index,2) -initnorms(index,1); -inittangs(index,2) inittangs(index,1)];
    revchangemat = [inittangs(index,1) initnorms(index,1); inittangs(index,2) initnorms(index,2)];

%loop over nearest neighbors and rewrite in terms of local coordinates    
     for i = 1:length(isites)                
         xtemp = isites(1,i);       
         ytemp = isites(2,i);
                
         isites(1,i) = changemat(1,1)*xtemp+changemat(1,2)*ytemp;
         isites(2,i) = changemat(2,1)*xtemp+changemat(2,2)*ytemp;
     end
            
%             if fig ==1
%                 figure
%                 title('Check Change of Basis')
%                 axis equal
%                 hold on
%                 plot(isites(1,:),isites(2,:), 'ko')
%                 pause
%             end
            
            xn = isites(1,:)';
            
            [X1,X2]=meshgrid(xn);
            
            A11 = phi1D(X1,X2,m);
            A12 = (xn(:,ones(1,p+1)).^x1D_exp(ones(1,n),:));
            
            Ax11 = phi1D_Dx(X1,X2,m);
            Ax12 = (xn(:,ones(1,p+1)).^x1D_exp_Dx(ones(1,n),:));
    
            B11 = isites(2,:);
            B12 = x1D_exp.*(xn(1,ones(1,p+1)).^x1D_exp);
                       
%             Bx11 = phi1Dx(X1(1,:),X2(1,:),m);
%             Bx12 = x1D_exp.*(xn(1,ones(1,p+1)).^x1D_exp_Dx);
            
            B = [B11 B12];
            
%            Bx = [Bx11 Bx12];
 
            A = [A11 A12;A12' zeros(p+1)];
            
          Ax = [Ax11 Ax12;Ax12' zeros(p+1)];

            S1D_local = B/A; %local interpolation weights
 
%            S1D(index,idx(index,:)) = S1D_local(1,1:n);        
            
            if fig ==2
                Xinterp = linspace(min(isites(1,:)),max(isites(2,:)),2*n)';
                [Xinterp1, Xinterp2] = meshgrid(isites(1,:),Xinterp);
                Ainterp11 = phi1D(Xinterp1,Xinterp2,m);
                Ainterp12 = (Xinterp(:,ones(1,p+1)).^x1D_exp(ones(1,2*n),:));
                Ainterp = [Ainterp11 Ainterp12];
                lambda = S1D_local';
                Yinterp = Ainterp*lambda;
                figure
                hold on
                axis equal
                title('local RBF interpolation')
                plot(isites(1,:), isites(2,:), 'ko')
                plot(Xinterp,Yinterp,'b-')
                pause
            end
            
            % Compute local norm direction
            localnorm = Ax*S1D_local';
            sclnorm = 1./sqrt(isites(1,1).^2+localnorm(1,1).^2);  %normalize scalar
            localnorm(1,1) = sclnorm.*localnorm(1,1); %scale vector
          
      if fig == 2
           figure
            hold on
            axis equal
            title('testing local normal direction')
            delta = 0.05;
            plot(isites(1,:), isites(2,:), 'ko')
            plot([isites(1,1) isites(2,1)],[isites(1,1) (isites(2,1)+delta*localnorm(1,1))],'r-')
  pause     
      end
      
            %revert basis of norms back to global coordinates
            calcnorms(index,1) = revchangemat(1,1)*isites(1,1) + revchangemat(1,2)*(isites(2,1)+localnorm(1,1));
            calcnorms(index,2) = revchangemat(2,1)*isites(1,1) + revchangemat(2,2)*(isites(2,1)+localnorm(1,1));
           
end
%% Normalizing the normal vectors at each node
% scal = sqrt(calcnorms(:,1).^2 + calcnorms(:,2).^2); %normalizing constant
% calcnorms(:,1) = calcnorms(:,1)./scal;  %normalize x-component
% calcnorms(:,2) = calcnorms(:,2)./scal;%normalize y-component
% calcnorms = calcnorms./scal; %normalize 

% plot initial calcualtion of normals 
% if fig == 1; 
% figure
%  title('Initial Calculated Normals Plot')
%  delta = 0.3;
%  hold on
%  axis equal
%  plot((X(1:N,1)),(X(1:N,2)),'ko')
%  for i = 1:N
%     plot([X(i,1) (X(i,1)+(delta*calcnorms(i,1)))],[X(i,2) (X(i,2)+(delta*calcnorms(i,2)))],'r')
%  end
% end

%% Correcting normal directions
%for i = 1:2N;
%    if  i<N || N+1<= i <= 2N-1;
for i = 1:N;
    if i < N;
      dir = sign((calcnorms(i,1)*calcnorms(i+1,1) + calcnorms(i,2)*calcnorms(i+1,2)));
        calcnorms(i+1,:) = dir.*calcnorms(i+1,:);        
    %elseif i == N || i == 2N;
    elseif i == N
        dir = sign((calcnorms(i,1)*calcnorms(1,1) + calcnorms(i,2)*calcnorms(1,2)));
        calcnorms(1,1) = dir.*calcnorms(1,1);  
        calcnorms(1,2) = dir.*calcnorms(1,2);
    end
end

calcnorms = -1.*calcnorms;

% % %Uncomment to plot direction normal check 
% if fig == 1
% figure
%  title('Corrected Calculated Normals Plot')
% delta = 0.3;
% hold on
% axis equal
% plot((X(1:N,1)),(X(1:N,2)),'ko')
% for i = 1:N
%    plot([X(i,1) (X(i,1)+(delta*calcnorms(i,1)))],[X(i,2) (X(i,2)+(delta*calcnorms(i,2)))],'r')
% end
% end

%% Compute trunorms
%evaluate parameterized curve from tensors at each gamma node 
t = linspace(0,2*pi,N+1)';
trunorms(:,1) = xtrunorms(t(1:N));
trunorms(:,2) = ytrunorms(t(1:N));
 scal = sqrt(trunorms(:,1).^2 + trunorms(:,2).^2); %normalizing constant
 trunorms(:,:) = trunorms(:,:)./scal; %normalize vector


% %plot the true normal directions
% if fig == 1
% figure
% title('True Normal Directions')
% delta = 0.3;
% hold on
% axis equal
% plot((X(1:N,1)),(X(1:N,2)),'ko')
% for i = 1:N
%    plot([X(i,1) (X(i,1)+(delta.*trunorms(i,1)))],[X(i,2) (X(i,2)+delta.*(trunorms(i,2)))],'g')
% end
% end

% %plot the true normal directions and calculated normals
% if fig == 1
% figure
% title('True Normal & Calculated Noraml Directions')
% delta = 0.3;
% hold on
% axis equal
% plot((X(1:N,1)),(X(1:N,2)),'ko')
% for i = 1:N
%    plot([X(i,1) (X(i,1)+(delta.*trunorms(i,1)))],[X(i,2) (X(i,2)+delta.*(trunorms(i,2)))],'g')
%    plot([X(i,1) (X(i,1)+(delta.*calcnorms(i,1)))],[X(i,2) (X(i,2)+delta.*(calcnorms(i,2)))],'r')
% end
% end


%% Compute relative normal error
errnormsx = norm(calcnorms(:,1) - trunorms(:,1),inf)/norm(calcnorms(:,1),inf);
errnormsy = norm(calcnorms(:,2) - trunorms(:,2),inf)/norm(calcnorms(:,2),inf);

errn1 = sqrt((trunorms(:,1) - calcnorms(:,1)).^2+(trunorms(:,2) - calcnorms(:,2)).^2)./sqrt((trunorms(:,1)).^2+(trunorms(:,2)).^2);
errn2 = sqrt((-trunorms(:,1) - calcnorms(:,1)).^2+(-trunorms(:,2) - calcnorms(:,2)).^2)./sqrt((trunorms(:,1)).^2+(trunorms(:,2)).^2);
errn = min(errn1,errn2);
errnorms = norm(errn,inf);

%Plot nomral error at each node
if fig == 1


figure
title('Log10 plot of L2 normal error at each node')
hold on 
axis equal
sz = 25;
scatter(X(1:N,1),X(1:N,2),sz,log10(errn),'filled')
colorbar

end

%% store nomral l2 errors
% errnorml2(nodes, neighnum-1) = log10(size(X,1));
% errnorml2(nodes +  size(nodedist,1),neighnum-1) = log10(errnorms);
errnorml2(nodes, polydeg- (minpolydeg - 1)) = size(X,1);
errnorml2(nodes +  size(nodedist,1),polydeg-(minpolydeg -1)) = errnorms;

    end
    
end

if fig == 3
%plot l2 error of normal 
figure
hold on
plot((log10(errnorml2(1:nodes,1))),(log10(errnorml2((nodes+1):2*nodes,1))),'b-o')
plot((log10(errnorml2(1:nodes,2))),(log10(errnorml2((nodes+1):2*nodes,2))),'k-o')
plot((log10(errnorml2(1:nodes,3))),(log10(errnorml2((nodes+1):2*nodes,3))),'g-o')
plot((log10(errnorml2(1:nodes,4))),(log10(errnorml2((nodes+1):2*nodes,4))),'r-o')
%plot((log10(errnorml2(1:nodes,5))),(log10(errnorml2((nodes+1):2*nodes,5))),'c-o')
%plot((log10(errnorml2(1:nodes,6))),(log10(errnorml2((nodes+1):2*nodes,6))),'m-o')
title('Normal line l2 error -- Increasing Neighbors and Polynomial Terms')
xlabel('Log10(N)')
ylabel('Log10(err)')
legend({'n = 6', 'n = 8', 'n = 10', 'n = 12'},'location','southwest')
end