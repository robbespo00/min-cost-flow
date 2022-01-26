% PRECONDITIONING: 
%   - d is a the vector corresponding to the diagonal of the matrix D;
%   - E is a sparse matrix that corresponds to the adjacency matrix of a
%       weighted directed graph;
%   - b1 is a vector that corresponds to the upper part of the known term;
%   - c1 is a vector that corresponds to the lower part of the known term;
%   - m is the number of iterations of the Arnoldi process;

function x = customGMRES(d, E, b1, c1, m)

    % we create the vector b that it's composed by b1 and c1
    b = [b1; c1];

    tic
    % The function "arnoldiqr" returns Q, R and Qn where:
    % - Q and R are the QR decomposition of the matrix Hn;
    % - Qn is the matrix from the Arnoldi process;
    [Q, R, Qn] = arnoldiqr(d, E, b, m);
    
    % We solve the least square problem using the the submatrices of Q
    % and R (it's cheaper)
    y = R(1:end-1,:) / Q(1, 1:end-1)'*norm(b);
    toc

    % We transform the solution of the smaller least square problem into
    % the solution of the bigger one inside the Krylov subspace
    x = Qn(:, 1:end-1)*y;
