% PRECONDITIONING: 
%   - diag is a the vector corresponding to the diagonal of the matrix D;
%   - E is a sparse matrix that corresponds to the incidence matrix of a
%       weighted directed graph;
%   - b1 is a vector that corresponds to the upper part of the known term;
%   - c1 is a vector that corresponds to the lower part of the known term;
%   - m is the number of iterations of the Arnoldi process;

function [x,q_nop, q_p,residual_nop, residual_p,t_nop, t_p2, t_px] = customGMRES(D, E, b1, c1, m, precond, a)

    edges = length(D);
    % we create the vector b that it's composed by b1 and c1
    b = [b1; c1];
   
    tic;
    [Q, R, Qn] = arnoldiqr(D, E, b, m);

    % We compute the solution of the system R_1*y=||b||*p_1^n (where p_1^n
    % is defined in the report)
    y_sol = R(1:end-1,:) \ (Q(1, 1:end-1)*norm(b))';

    % We compute the solution of the main problem
    x = Qn(:, 1:end-1)*y_sol;

    Ax = [D.*x(1:length(D)) + E'*x(length(D)+1:end); E*x(1:length(D))];
    
    residual_nop = norm(Ax-b)/norm(b);
        
    t_nop = toc;
    q_nop=abs(Q(1,m+1));
    
    
    if precond
        D3 = D.^3; 
        b1 = D.*b1 + E'*c1;
        c1 = -E*((E'*c1)./D)+ a*c1;
        b = [b1; c1]; % P^T*b
        tic;
        [Q, R, Qn] = precondarnoldi(D, E, b, m, a);
        y_sol = R(1:end-1,:) \ (Q(1, 1:end-1)*norm(b))';
        y = Qn(:, 1:end-1)*y_sol; % sol. of P^T*A*P*y=P^T*b_tilde 
        t_p = toc;
        
        y1 = y(1:edges);
        y2 = y(edges+1:end);
        
        
        tic;
        z1 = (D3.*y1) + E'*(E*(D.*y1)) + D.*(E'*(E*y1)) -D.*(E'*(E*((E'*y2)./D)));
        z2 = -E*((E'*(E*(D.*y1)))./D);
        z=[z1;z2];
        z = [z(1:edges) + a*(D.*(E'*y2)); z(edges+1:end)+ a*(E*(D.*y1))]; %P^T*A*P*y
        residual_p = norm(z-b)/norm(b);
        t_p2= toc;
        t_p2 = t_p2+t_p;
        
        q_p=abs(Q(1,m+1)); % lower bound of Preconditioned Transformed
        
        % Preconditioned Original x=P*y
        tic;
        x1 = D.*y1;
        x2 = E*y1-E*((E'*y2)./D) + a*y2;
        x = [x1; x2];
        t_px = toc;
        t_px = t_px + t_p;
    end
        

end