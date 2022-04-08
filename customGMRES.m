% PRECONDITIONING: 
%   - diag is a the vector corresponding to the diagonal of the matrix D;
%   - E is a sparse matrix that corresponds to the incidence matrix of a
%       weighted directed graph;
%   - b1 is a vector that corresponds to the upper part of the known term;
%   - c1 is a vector that corresponds to the lower part of the known term;
%   - m is the number of iterations of the Arnoldi process;

function [x,q_nop, q_p,residual_nop, residual_p,t_nop, t_p2, t_px] = customGMRES(D, E, b1, c1, m, precond, a)

    nodes = size(E,1);
    edges = length(D);
    % we create the vector b that it's composed by b1 and c1
    b = [b1; c1];
    D3 = D.^3; 
    tic;
    %if precond
        [Q, R, Qn] = arnoldiqr(D, E, b, m);
        q_nop=abs(Q(1,m+1));
    
    % We compute the solution of the system R_1*y=||b||*p_1^n (where p_1^n
    % is defined in the report)
    y = R(1:end-1,:) \ (Q(1, 1:end-1)*norm(b))';
    

    % We compute the solution of the main problem
    x = Qn(:, 1:end-1)*y;
        
    
    
        Ax = [D.*x(1:length(D)) + E'*x(length(D)+1:end); E*x(1:length(D))];
        residual_nop = norm(Ax-b)/norm(b);
        
    t_nop = toc;
        
    
        b1 = D.*b1 + E'*c1;
        c1 = -E*((E'*c1)./D)+ a*c1;
        b = [b1; c1];
       
    tic;
        [Q, R, Qn] = precondarnoldi(D, E, b, m, a);
        
    %else
        
    %end

    % The function "arnoldiqr" returns Q, R and Qn where:
    % - Q and R are the QR decomposition of the matrix Hn;
    % - Qn is the matrix from the Arnoldi process;
    
       
    q_p=abs(Q(1,m+1));
    
    % We compute the solution of the system R_1*y=||b||*p_1^n (where p_1^n
    % is defined in the report)
    y = R(1:end-1,:) \ (Q(1, 1:end-1)*norm(b))';
    

    % We compute the solution of the main problem
    x = Qn(:, 1:end-1)*y;
    qk1 = x(1:edges);
    qk2 = x(edges+1:end);
    
     t_p = toc;
    
    if precond

        tic;
        z1 = (D3.*qk1) + E'*(E*(D.*qk1)) + D.*(E'*(E*qk1)) -D.*(E'*(E*((E'*qk2)./D)));
        z2 = -E*((E'*(E*(D.*qk1)))./D);
        z=[z1;z2];
        z = [z(1:edges) + a*(D.*(E'*qk2)); z(edges+1:end)+ a*(E*(D.*qk1))];
        

       residual_p = norm(z-b)/norm(b)
       t_p2= toc;
       t_p2 = t_p2+t_p;

       tic;
       x1 = D.*x(1:edges);
       x2 = E*x(1:edges)-E*((E'*x(edges+1:end))./D) + a*x(edges+1:end);
       x = [x1; x2];
       t_px = toc;
       t_px = t_px + t_p;
%         D2 = D.^2;
%         x1 = x(1:edges);
%         x2 = x(edges+1:end);
% 
%         PtAx = [D2.*x1 + E'*(E*x1) + D.*(E'*x2); -E*((E'*(E*x1))./D) + a*(E*x1)];
% 
%         residual = norm(PtAx-b)/norm(b)
        
    end