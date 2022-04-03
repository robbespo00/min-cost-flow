% PRECONDITIONING: 
%   - diag is a the vector corresponding to the diagonal of the matrix D;
%   - E is a sparse matrix that corresponds to the incidence matrix of a
%       weighted directed graph;
%   - b1 is a vector that corresponds to the upper part of the known term;
%   - c1 is a vector that corresponds to the lower part of the known term;
%   - m is the number of iterations of the Arnoldi process;

function [x,q] = customGMRES(D, E, b1, c1, m, precond, a)
    
    nodes = size(E,1);
    edges = length(D);
    % we create the vector b that it's composed by b1 and c1
    b = [b1; c1];
    
    if precond
        b1 = D.*b1 + E'*c1;
       
        c1 = -E*((E'*c1)./D)+ a*c1;
        
        b = [b1; c1];
       
        [Q, R, Qn] = precondarnoldi(D, E, b, m, a);
        
    else
        [Q, R, Qn] = arnoldiqr(D, E, b, m);
    end

    % The function "arnoldiqr" returns Q, R and Qn where:
    % - Q and R are the QR decomposition of the matrix Hn;
    % - Qn is the matrix from the Arnoldi process;
    
       
    q=abs(Q(1,m+1));
    
    % We compute the solution of the system R_1*y=||b||*p_1^n (where p_1^n
    % is defined in the report)
    y = R(1:end-1,:) \ (Q(1, 1:end-1)*norm(b))';
    

    % We compute the solution of the main problem
    x = Qn(:, 1:end-1)*y;
    
    S = -E*(E'./D);
    P = [diag(D) zeros(edges, nodes); E S];

    A = [diag(D) E'; E zeros(nodes, nodes)];
    B = P'*A*P;
    save('B','B');
    save('b1','b');
    By = B*x;
    
    residual = norm(By-b)/norm(b)
    
    
    if precond
        x1 = D.*x(1:edges);
        x2 = E*x(1:edges)-E*((E'*x(edges+1:end)./D)) + a*x(edges+1:end);
        x = [x1; x2];
    end

