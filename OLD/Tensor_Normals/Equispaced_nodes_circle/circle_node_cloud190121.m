%Evenly spaced node distribution on unit circle and with cloud.

close all
clear all

%---Inputs---
xran = [-1.5 1.5]; yran = [-1.5 1.5]; % x and y range
maxgN = 2080;  %maximum number of cloud nodes

%---Parameterized Unit Circle---
xfunc = @(th) cos(th);
yfunc = @(th) sin(th);


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









