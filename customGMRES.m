% PRECONDITIONING: 
%   - diag is a the vector corresponding to the diagonal of the matrix D;
%   - E is a sparse matrix that corresponds to the incidence matrix of a
%       weighted directed graph;
%   - b1 is a vector that corresponds to the upper part of the known term;
%   - c1 is a vector that corresponds to the lower part of the known term;
%   - m is the number of iterations of the Arnoldi process;

function x = customGMRES(diag, E, b1, c1, m)

    % we create the vector b that it's composed by b1 and c1
    b = [b1; c1];

    
    % The function "arnoldiqr" returns Q, R and Qn where:
    % - Q and R are the QR decomposition of the matrix Hn;
    % - Qn is the matrix from the Arnoldi process;
    [Q, R, Qn] = arnoldiqr(diag, E, b, m);

    % We compute the solution of the system R_1*y=||b||*p_1^n (where p_1^n
    % is defined in the report)
    y = R(1:end-1,:) \ (Q(1, 1:end-1)*norm(b))';
    

    % We compute the solution of the main problem
    x = Qn(:, 1:end-1)*y;

