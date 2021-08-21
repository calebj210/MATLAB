%%Generates Hexagonal Grid
%%N: Total number of nodes approximately
%%bounding_box: [left right down up] 
function [X Y] = hex(N,bounding_box)
Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(0:1:round(sqrt(N))-1);
n = size(X,1);
% X = Rad3Over2 * X;
X = X/max(X(:));
Y = Y + repmat([0 0.5],[n,n/2]);
Y = Y/max(Y(:));
Lx=bounding_box(2)-bounding_box(1);
Ly=bounding_box(4)-bounding_box(3);
X=X(:)*Lx+bounding_box(1);
Y=Y(:)*Ly+bounding_box(3);