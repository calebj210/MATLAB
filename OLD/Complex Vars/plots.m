% Test code for displaying an analytic function; here exp(z)
clear; close all;
bx = [-4,4,-10,10]; % Domain (box) surrounding origin to be displayed
nx = 41; % Number of nodes along real axis for re & im plots
ny = 41; % Number of nodes along imag axis for re & im plots
bounds = [ -60,60; -60,60; 0,60]; % Lower and upper in the three displays
f = @(z) z.^2; % Function to display
display_function(f,bx,nx,ny,bounds); % Create the three subplots

