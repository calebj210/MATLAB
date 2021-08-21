a = -2;
b = 5;
n = 10;
f=-(linspace(a,b,n)).^3.+3;
trapz(linspace(a,b,n),f)