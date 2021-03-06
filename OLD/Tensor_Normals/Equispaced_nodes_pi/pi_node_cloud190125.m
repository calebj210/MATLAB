%Evenly spaced node distribution on pi surface and with cloud.

close all
clear all

%---Inputs---
xran = [-1.5 1.5]; yran = [-1.5 1.5]; % x and y range
maxgN = 10080;  %maximum number of cloud nodes

% %---Parameterized Unit Circle---
% xfunc = @(th) cos(th);
% yfunc = @(th) sin(th);


%---Parameterized Pi Curve ---
xfunc=@(t) ((17/31).*sin((235/57)+(-32).*t)+(19/17).*sin((192/55)+(-30).*t)+(47/32).*sin((69/25)+(-29).*t)+(35/26).*sin((75/34)+(-27).*t)+(6/31).*sin((23/10)+(-26).*t)+(35/43).*sin((10/33)+(-25).*t)+(126/43).*sin((421/158)+(-24).*t)+(143/57).*sin((35/22)+(-22).*t)+(106/27).*sin((84/29)+(-21).*t)+(88/25).*sin((23/27)+(-20).*t)+(74/27).*sin((53/22)+(-19).*t)+(44/53).*sin((117/25)+(-18).*t)+(126/25).*sin((88/49)+(-17).*t)+(79/11).*sin((43/26)+(-16).*t)+(43/12).*sin((41/17)+(-15).*t)+(47/27).*sin((244/81)+(-14).*t)+(8/5).*sin((79/19)+(-13).*t)+(373/46).*sin((109/38)+(-12).*t)+(1200/31).*sin((133/74)+(-11).*t)+(67/24).*sin((157/61)+(-10).*t)+(583/28).*sin((13/8)+(-8).*t)+(772/35).*sin((59/16)+(-7).*t)+(3705/46).*sin((117/50)+(-6).*t)+(862/13).*sin((19/8)+(-5).*t)+(6555/34).*sin((157/78)+(-3).*t)+(6949/13).*sin((83/27)+(-1).*t)+(-6805/54).*sin((1/145)+2.*t)+(-5207/37).*sin((49/74)+4.*t)+(-1811/58).*sin((55/43)+9.*t)+(-63/20).*sin((2/23)+23.*t)+(-266/177).*sin((13/18)+28.*t)+(-2/21).*sin((7/16)+31.*t))/1000;
yfunc=@(t) ((70/37).*sin((65/32)+(-32).*t)+(11/12).*sin((98/41)+(-31).*t)+(26/29).*sin((35/12)+(-30).*t)+(54/41).*sin((18/7)+(-29).*t)+(177/71).*sin((51/19)+(-27).*t)+(52/33).*sin((133/52)+(-26).*t)+(59/34).*sin((125/33)+(-26).*t)+(98/29).*sin((18/11)+(-25).*t)+(302/75).*sin((59/22)+(-24).*t)+(104/9).*sin((118/45)+(-22).*t)+(52/33).*sin((133/52)+(-21).*t)+(37/45).*sin((61/14)+(-20).*t)+(143/46).*sin((144/41)+(-19).*t)+(254/47).*sin((19/52)+(-18).*t)+(246/35).*sin((92/25)+(-17).*t)+(722/111).*sin((176/67)+(-16).*t)+(136/23).*sin((3/19)+(-15).*t)+(273/25).*sin((32/21)+(-13).*t)+(229/33).*sin((117/28)+(-12).*t)+(19/4).*sin((43/11)+(-11).*t)+(135/8).*sin((23/10)+(-10).*t)+(205/6).*sin((33/23)+(-8).*t)+(679/45).*sin((55/12)+(-7).*t)+(101/8).*sin((11/12)+(-6).*t)+(2760/59).*sin((40/11)+(-5).*t)+(1207/18).*sin((21/23)+(-4).*t)+(8566/27).*sin((39/28)+(-3).*t)+(12334/29).*sin((47/37)+(-2).*t)+(15410/39).*sin((185/41)+(-1).*t)+(-596/17).*sin((3/26)+9.*t)+(-247/28).*sin((25/21)+14.*t)+(-458/131).*sin((21/37)+23.*t)+(-41/36).*sin((7/8)+28.*t))/1000;
xpfunc =@(t) (-544/31).*cos((235/57)+(-32).*t)+(-570/17).*cos((192/55)+(-30).*t)+(-1363/32).*cos((69/25)+(-29).*t)+(-945/26).*cos((75/34)+(-27).*t)+(-156/31).*cos((23/10)+(-26).*t)+(-875/43).*cos((10/33)+(-25).*t)+(-3024/43).*cos((421/158)+(-24).*t)+(-3146/57).*cos((35/22)+(-22).*t)+(-742/9).*cos((84/29)+(-21).*t)+(-352/5).*cos((23/27)+(-20).*t)+(-1406/27).*cos((53/22)+(-19).*t)+(-792/53).*cos((117/25)+(-18).*t)+(-2142/25).*cos((88/49)+(-17).*t)+(-1264/11).*cos((43/26)+(-16).*t)+(-215/4).*cos((41/17)+(-15).*t)+(-658/27).*cos((244/81)+(-14).*t)+(-104/5).*cos((79/19)+(-13).*t)+(-2238/23).*cos((109/38)+(-12).*t)+(-13200/31).*cos((133/74)+(-11).*t)+(-335/12).*cos((157/61)+(-10).*t)+(-1166/7).*cos((13/8)+(-8).*t)+(-772/5).*cos((59/16)+(-7).*t)+(-11115/23).*cos((117/50)+(-6).*t)+(-4310/13).*cos((19/8)+(-5).*t)+(-19665/34).*cos((157/78)+(-3).*t)+(-6949/13).*cos((83/27)+(-1).*t)+(-6805/27).*cos((1/145)+2.*t)+(-20828/37).*cos((49/74)+4.*t)+(-16299/58).*cos((55/43)+9.*t)+(-1449/20).*cos((2/23)+23.*t)+(-7448/177).*cos((13/18)+28.*t)+(-62/21).*cos((7/16)+31.*t);
ypfunc =@(t) (-2240/37).*cos((65/32)+(-32).*t)+(-341/12).*cos((98/41)+(-31).*t)+(-780/29).*cos((35/12)+(-30).*t)+(-1566/41).*cos((18/7)+(-29).*t)+(-4779/71).*cos((51/19)+(-27).*t)+(-1352/33).*cos((133/52)+(-26).*t)+(-767/17).*cos((125/33)+(-26).*t)+(-2450/29).*cos((18/11)+(-25).*t)+(-2416/25).*cos((59/22)+(-24).*t)+(-2288/9).*cos((118/45)+(-22).*t)+(-364/11).*cos((133/52)+(-21).*t)+(-148/9).*cos((61/14)+(-20).*t)+(-2717/46).*cos((144/41)+(-19).*t)+(-4572/47).*cos((19/52)+(-18).*t)+(-4182/35).*cos((92/25)+(-17).*t)+(-11552/111).*cos((176/67)+(-16).*t)+(-2040/23).*cos((3/19)+(-15).*t)+(-3549/25).*cos((32/21)+(-13).*t)+(-916/11).*cos((117/28)+(-12).*t)+(-209/4).*cos((43/11)+(-11).*t)+(-675/4).*cos((23/10)+(-10).*t)+(-820/3).*cos((33/23)+(-8).*t)+(-4753/45).*cos((55/12)+(-7).*t)+(-303/4).*cos((11/12)+(-6).*t)+(-13800/59).*cos((40/11)+(-5).*t)+(-2414/9).*cos((21/23)+(-4).*t)+(-8566/9).*cos((39/28)+(-3).*t)+(-24668/29).*cos((47/37)+(-2).*t)+(-15410/39).*cos((185/41)+(-1).*t)+(-5364/17).*cos((3/26)+9.*t)+(-247/2).*cos((25/21)+14.*t)+(-10534/131).*cos((21/37)+23.*t)+(-287/9).*cos((7/8)+28.*t);


%---Node Distribution---
box = [xran yran]; % bounding box

%generate background nodes
[X,Y] = hex(maxgN,box); % hexagonal grid for background nodes


figure
hold on
plot(X,Y,'x')

%compute average distance of nodes in node cloud
[idx,dist] = knnsearch([X,Y],[X,Y],'K',7);

% %uncomment to plot nearest neighbors
% N = length(x);
% figure
% axis equal
% for i = 1:N
% hold on
% plot(x,y,'w*')
% plot(x,y,'ko')
% plot((x(N+1:length(x(:,1)),1)),(y(N+1:length(y(:,1)),1)),'b*')
% plot((x(idx(i,:),1)),(y(idx(i,:),1)),'r*')
% plot((x(i,1)),(y(i,1)),'g*')
% if i >1
% for j = 1:i-1
%     plot((x(j,1)),(y(j,1)),'mx')
% end
% end
% pause
% end


hdist = mean(dist);
h = (hdist(:,2)+hdist(:,3))/2;

cN = round((2*pi)/h);  %number of nodes on curve

t = linspace(0,2*pi,cN);

%storage for nodes on curve

gamma = zeros(cN,2);

gamma(:,1) = xfunc(t);
gamma(:,2) = yfunc(t);

figure
axis equal
hold on
plot(X,Y,'rx')
plot(gamma(:,1),gamma(:,2),'bo')

clear idx dist

[idx,dist]=knnsearch([X,Y],gamma,'K',50);




for i = 1:size(dist,1)
   for j = 1:size(dist,2)
      if dist(i,j) < 0.75*h
         X(idx(i,j),:) = 90;
         Y(idx(i,j),:) = 90;
         idx(i,j) = NaN;
         dist(i,j) = NaN;
      end
   end
end

kf = size(X,1);
for k = 1:(kf)
    
    if X(kf+1-k,1) == 90
        X(kf+1-k,:) = [];
        Y(kf+1-k,:) = [];
    end
% if X(k,1) == 90
%     X(k,:) =[];
%     Y(k,:) = [];
% end
end



figure
axis equal
hold on
plot(X,Y,'rx')
plot(gamma(:,1),gamma(:,2),'bo')









