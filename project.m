function [x,y,Hn,b_tilde] = project(D, E, b1, c1, m)

d_dim=size(D,1);
A=[D E'; E zeros(d_dim, d_dim)];
b=[b1; c1];
[Q,R, Qn,Hn]=arnoldiqr(A,b,m);
b_tilde=norm(b)*[1; zeros(length(Q)-1,1)];
%y=lsqr(R(1:end-1,:), Q(:,1:end-1)'*b_tilde,1e-15,600);
y=lsqr(Hn, b_tilde, 1e-15, 1000);

tot=Qn*Hn-A*Qn(:, 1:end-1);
norm(tot)
x=Qn(:,1:end-1)*y;
