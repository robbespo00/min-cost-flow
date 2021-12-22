% PRECONDITION: 
%   - d is a vector that compose the diagonal matrix D;
%   - E is a sparse matrix that corresponds to the adjacency matrix of a
%       weighted directed graph;
%   - b1 is a vector that corresponds to the first part of known terms;
%   - c1 is a vector that corresponds to the second part of known terms;
%   - m is an integer that corresponds to the number of Arnoldi iterations;
function [x] = project(d, E, b1, c1, m)

% we create the vector b that it's composed by b1 and c1
b = [b1; c1];

tic
% The function "arnoldiqr" returns Q, R and Qn where:
% - Q and R are the QR factorization of the matrix Hn;
% - Qn
[Q, R, Qn]=arnoldiqr(d, E, b, m);
b_tilde=norm(b)*[1; zeros(length(Q)-1,1)];
% We solve the problem of the least square using the the submatrices of Q
% and R such that the cost is cheaper
y=lsqr(R(1:end-1,:), Q(:,1:end-1)'*b_tilde,1e-15,5000);
toc
% We transform the solution of the cheaper least square into the original
% space
x=Qn(:,1:end-1)*y;