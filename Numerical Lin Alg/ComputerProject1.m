format short
format loose
clear
clc

%% Problem a %%
disp('%%%%% Probelm a %%%%%')
A1=[2,1,1,0;...
	   4,3,3,1;...
	   8,7,9,5;...
	   6,7,9,8];
b1=[-1;3;2;1];

GEPP(A1,b1);

%% Problem b %%
disp('%%%%% Probelm b %%%%%')
A2=[5,-7,-4,9;...
	   6,-8,-7,5;...
	   4,-4,-9,-9;...
	  -9,11,16,7];
b2=[2;-1;0;1];

GEPP(A2,b2);


%% Problem c %%
disp('%%%%% Probelm c %%%%%')
A3i=[1,1./2,1./3,1./4;...
	    1./2,1./3,1./4,1./5;...
		1./3,1./4,1./5,1./6;...
		1./4,1./5,1./6,1./7];
A3ii=[1.00,0.50,0.33,0.25;...
	     0.50,0.33,0.25,0.20;
		 0.33,0.25,0.20,0.16;
		 0.25,0.20,0.16,0.14];
b3=[1;0;0;0];

% c.i
disp('%% c.i %%')
GEPP(A3i,b3);

% c.ii
disp('%% c.ii %%')
GEPP(A3ii,b3);