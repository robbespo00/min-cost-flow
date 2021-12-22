function [Q, R, Qn] = arnoldiqr(d, E, q1, m)

    n = 2*length(d);
    q1 = q1/norm(q1);
    Qn = zeros(n, m+1); Qn(:, 1) = q1;
    Hn = zeros(m+1, m);
    Q = zeros(m+1, m+1);
    R = zeros(m+1, m); 
    
    for k = 1:m
        z = [ (d.*Qn(1:n/2, k)) + (E'*Qn((n/2)+1:n, k)); E*Qn(1:n/2, k)]; %we want to avoid full computation since the matrices are sparse
        
        for i = max(k-1, 1):k
            Hn(i, k) = Qn(:, i)'*z;
            z = z - Hn(i, k)*Qn(:, i);
        end
        
        Hn(k+1, k) = norm(z);
        
        if Hn(k+1, k) < 1e-10   %breakdown% 
            display("BREAKDOWN HAPPENED!");
            return; 
        end 
        
        Qn(:, k+1) = z/Hn(k+1, k);
        
        if k == 1
            x = Hn(1:2, 1);
            y = norm(x)*[1; 0];
            v = x - y;
            v_norm = 2/(norm(v)^2);
            vHn = v'*Hn(1:2, 1);
            T = v.*vHn;
            R(1:2, 1) = Hn(1:2, 1) - v_norm*T;
            Q(1:2, 1:2) = (eye(2)-v_norm*v.*v')';
        else 
            %Building R tilde%
            c = Hn(k-1, k)*Q(k-1, :)' + Hn(k, k)*Q(k, :)';
            x = [c(k); Hn(k+1, k)];
            x_norm = norm(x);
            y = x_norm*[1; 0];   
            v = x - y;
            v_norm = norm(v)^2;
            R(1:k, k) = [c(1:k-1); x_norm];
            
            %Building Q tilde%
            Q(k+1,k+1)=1; %build Q' segnato%
            u = Q(1:k-1, k);
            a = Q(k, k);
            p = [u zeros(k-1, 1); a 0; 0 1]; %last two rows of Q segnato transponed %
            w = [v(1)^2*u v(1)*v(2)*u; v(1)^2*a v(1)*v(2)*a; v(1)*v(2) v(2)^2];
            Q(1:k+1, k:k+1) = p-(2/(v_norm))*w;
        end
    end
    
    



