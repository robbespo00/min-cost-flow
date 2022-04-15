function [Q, R, Qn] = arnoldiqr(D, E, q1, m, a, precond) % initialization
    
    d = length(D);
    tau = size(E,1);
    n = d + tau;

    if m <= 0
        display("Error. m must be > 0");
        return;
    elseif m > n
        m = n;
        display("You don't need more than n iterations. m has been set to n");
    end

    q1 = q1/norm(q1);
    Qn = zeros(n, m+1); 
    Qn(:, 1) = q1;
    Hn1 = zeros(m,1); %diagonal of H_n
    Hn2 = zeros(m,1); %subdiagonal of H_n
    Q = zeros(m+1, m+1);
    R = zeros(m+1, m);

    if precond
        D3 = D.^3; 
    end

    for k = 1:m
        
        qk1 = Qn(1:d, k);
        qk2 = Qn(d+1:n, k);

        if precond %generating the new vector of the Krylov Subspace (precond or not)
            z1 = (D3.*qk1) + E'*(E*(D.*qk1)) + D.*(E'*(E*qk1)) -D.*(E'*(E*((E'*qk2)./D)));
            z2 = -E*((E'*(E*(D.*qk1)))./D);
            z = [z1; z2];
    
            z = [z(1:d) + a*(D.*(E'*qk2)); z(d+1:end)+ a*(E*(D.*qk1))];
        else
            z = [ (D.*qk1) + (E'*qk2); E*qk1]; 
        end

        %COMPUTING NEW COMPONENTS OF Hn
        if k==1
            Hn1(1) = Qn(:,1)'*z;
            z = z - Hn1(1)*Qn(:, 1);
            Hn2(1) = norm(z);
        else
            z = z - Hn2(k-1)*Qn(:, k-1);
            Hn1(k)= Qn(:, k)'*z;
            z = z - Hn1(k)*Qn(:, k);
            Hn2(k) = norm(z);
        end

        if abs(Hn2(k)) < 1e-10 % check for breakdown
            disp("BREAKDOWN HAPPENED!");
            return;
        end

        Qn(:, k+1) = z/Hn2(k); %new element of the orthonormal basis of Krylov Subspace

        %QR DECOMPOSITION OF Hn (1st AND OTHER ITERATIONS)
        if k == 1
            x = [Hn1(1); Hn2(1)];
            y = norm(x)*[1; 0];
            v = x - y;
            v_norm = 2/(norm(v)^2);
            vHn = v'*x;
            T = v.*vHn;
            R(1:2, 1) = x - v_norm*T;
            Q(1:2, 1:2) = (eye(2)-v_norm*v.*v')';
        else
            %Building R tilde%
            c = Hn2(k-1)*Q(k-1, 1:k)' + Hn1(k)*Q(k, 1:k)';
            x = [c(k); Hn2(k)];
            x_norm = norm(x);
            y = x_norm*[1; 0];
            v = x - y;
            v_norm = norm(v)^2;
            R(1:k, k) = [c(1:k-1); x_norm];

            %Building Q tilde%
            Q(k+1,k+1)=1; %build \underline{Q}^T
            u = Q(1:k-1, k);
            a2 = Q(k, k);
            p = [u zeros(k-1, 1); a2 0; 0 1]; %last two rows of \underline{Q}^T
            gamma = [v(1)^2*u v(1)*v(2)*u; v(1)^2*a2 v(1)*v(2)*a2; v(1)*v(2) v(2)^2];
            Q(1:k+1, k:k+1) = p-(2/(v_norm))*gamma;
        end
    end
end
