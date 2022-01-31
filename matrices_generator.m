function [A, E, D, b, c] = matrices_generator(n, e)

    E=graph_generating(n,e);
    D = 5*rand(e,1);
    b= 10*(rand(e,1)-(ones(e,1))/2);
    c= 10*(rand(n,1)-(ones(n,1))/2);
    A=[D.*eye(e) full(E)'; full(E) zeros(n,n)];