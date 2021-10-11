function flipped_day4

% this code compares the performance of Hermite and Lagrange
% interpolation for two functions

% f = @(x) exp(x).*sin(x);
% fp = @(x) exp(x).*cos(x)+exp(x).*sin(x);

f = @(x) exp(-16*x.^2);
fp = @(x) -32*x.*exp(-16*x.^2);

xeval = [-1:0.01:1];

for N =2:2:20

x = linspace(-1,1,N);
y = f(x);
yp = fp(x);



yL = lagrange(xeval,x,y);

yH = hermite(xeval,x,y,yp);

figure(1)
plot(xeval,yL,xeval,yH,xeval, f(xeval),'LineWidth',3)
legend('Lagrange','Hermite','f(x)')
title('Plots of approximations')

figure(2)
semilogy(xeval,abs(yL-f(xeval)),xeval,abs(yH-f(xeval)),'LineWidth',3)
legend('Lagrange','Hermite')
title('Absolute error')

    
pause
end


return

function y = hermite(x,pointx,pointy, pointyp)

% hermite interpolation
n=size(pointx,2);
L=ones(n,size(x,2));
if (size(pointx,2)~=size(pointy,2))
    fprintf(1,'\nERROR!\nPOINTX and POINTY must have the same number of elements\n');
    y=NaN;
else
    % create matrix with Lagrange intepolation
    for i=1:n
        for j=1:n
            if (i~=j)
                L(i,:)=L(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
            end
        end
    end
    
    %   create the constants L_j'(x_j) needed for Q_j
    ll = zeros(n,1);
    
    for j = 1:n
        for l = 1:n
            if (j~=l)
                ll(j) = ll(j)+1/(pointx(j)-pointx(l));
            end
        end
    end
    
    % Now evaluate R_j and Q_j
    
    Q = zeros(size(L));
    R = zeros(size(L));
    
    for j = 1:n
        Q(j,:) = (1-2*(x-pointx(j))*ll(j)).*(L(j,:)).^2;
        R(j,:) = (x-pointx(j)).*(L(j,:)).^2;
    end
    
    % finally evaluate the polynomial
    
    y = zeros(size(x));
    
    for j = 1:n
        y = y + pointy(j)*Q(j,:)+pointyp(j)*R(j,:);
    end
    
end

return


function y=lagrange(x,pointx,pointy)
%
%LAGRANGE   approx a point-defined function using the Lagrange polynomial interpolation
%
%      LAGRANGE(X,POINTX,POINTY) approx the function definited by the points:
%      P1=(POINTX(1),POINTY(1)), P2=(POINTX(2),POINTY(2)), ..., PN(POINTX(N),POINTY(N))
%      and calculate it in each elements of X
%
%      If POINTX and POINTY have different number of elements the function will return the NaN value
%
%      function wrote by: Calzino
%      7-oct-2001
%
n=size(pointx,2);
L=ones(n,size(x,2));
if (size(pointx,2)~=size(pointy,2))
    fprintf(1,'\nERROR!\nPOINTX and POINTY must have the same number of elements\n');
    y=NaN;
else
    for i=1:n
        for j=1:n
            if (i~=j)
                L(i,:)=L(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
            end
        end
    end
    y=0;
    for i=1:n
        y=y+pointy(i)*L(i,:);
    end
end

return