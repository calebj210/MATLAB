%{
% Solve generalized eigenvalue problem.
% Authors: Alan Bouwman, Caleb Jacobs
%}

% Clean workspace
clear;
close all;

% Load dense matrices
A1 = load('Test_Matrices/A1.mat');
A2 = load('Test_Matrices/A2.mat');

A1 = A1.Expression1;
A2 = A2.Expression1;

% Load sparse matrices using our mtx reader defined below
B1 = readMTX(load('Test_Matrices/B1.mtx'));
B2 = readMTX(load('Test_Matrices/B2.mtx'));

% Solve generalized eigenvalue problems
[VA, DA] = eigs(A1, A2, length(A1)); % Eigen system of A1*x = lambda*A2*x
[VB, DB] = eigs(B1, B2, length(B1)); % Eigen system of B1*x = lambda*B2*x

% Create vectors of eigenvalues
lamA = diag(DA);
lamB = diag(DB);

% Function for reading mtx files into a sparse matrix
function A = readMTX(Atmp)
    % Pass row, column, and value pairs to the sparse function
    A = sparse(Atmp(2:end,1), ...
               Atmp(2:end,2), ...
               Atmp(2:end,3), ...
               Atmp(1,1), ...
               Atmp(1,2));
end
