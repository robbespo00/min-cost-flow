function init_customGMRES(filename, mode, generate, distribution, precond, a)
    
    prompt="Insert the number of nodes of the graph: ";
    nodi = input(prompt);
    if ~isnumeric(nodi)
        ME = MException('%s is not a number',nodi);
        throw(ME)
    end
    [E, b, c] = netgenreader(filename, nodi);
    e = size(E,2);
    
    % given b and c we build b tilde
    b=-b;
    b_tilde = [b; c];
    

    % compute the norm of b tilde
    b_norm = norm(b_tilde);
    n = nodi + e;
    prompt="Insert delta value. Insert 0 for the default value: ";
    in = input(prompt);
    if in == 0
        delta = 10^(-10);
    else 
        if isnumeric(in) && in > 0
            delta = in;
        else
            ME = MException('%s is not a number or is negative',in);
            throw(ME)
        end
    end

    if generate
        rng(46273);
        switch distribution
        case 'Gamma'
            D = random('Gamma', 5,1, e, 1);
            save('D_gamma', 'D');
        case 'Beta75'
            D = random('Beta', 0.75, 0.75, e, 1);
            save('D_beta75', 'D');
        case 'Beta44'
            D = random('Beta', 4, 4, e, 1);
            save('D_beta44', 'D');
        case 'Chi'
            D = random('chi2', 4, e, 1);
            save('D_chi2', 'D');
        case 'MixBin'
             D=random('Binomial', 20, 0.7, e,1)+9;
             perm = randperm(e, round(e*0.7));
             D(perm) = random('Binomial', 20, 0.7, round(e*0.7),1);
             save('D_mixbin', 'D');
        case 'Uniform'
            D = rand(e,1);
            save('D_uniform', 'D');
        case 'Ill'
             D=rand(e,1);
             perm = randperm(e, round(e*0.4));
             perm = randperm(e, round(e*0.4));
             D(perm)= D(perm)*10^-5;
             save('D_ill', 'D');
        end
    else
        load(distribution)
    end
    

    switch mode
        case 'minres'
            % Building A in sparse format
            A = sparse([1:e], [1:e], D, nodi+e, nodi+e);
            A(e+1:end,1:e) = E;
            A(1:e, e+1:end) = E';
            % Create the file where the results will be saved
            fileID = fopen(strcat(erase(filename, ".txt"),"_result.txt"), 'w');
            jump = round(sqrt(n)/5);
        case 'precond'
            jump = round(sqrt(n)/5);
        otherwise
            jump=1;
    end

    
    
    num_iterations = [50:jump:round(n/25)];

    
    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES with different number of iterations 
    t_nop = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES with different number of iterations
    residual_nop = zeros(length(num_iterations), 1);
    
    if strcmp(mode, 'rate')
        rate = zeros(length(num_iterations)-1, 1);
    end
    
    q_nop = zeros(length(num_iterations)-1, 1);
    
    if precond
        residual_p = zeros(length(num_iterations),1);
        flag=true;
        residual_px = zeros(length(num_iterations),1);

        t_p2 = zeros(length(num_iterations),1);
        t_px = zeros(length(num_iterations),1);
        q_p = zeros(length(num_iterations),1);
    end
    
    
    for i = 1:length(num_iterations)
        fprintf('Iteration number %d (m=%d)\n', i, num_iterations(i));
       
        
        if precond 
            [x_custom(:, i), q_nop(i),q_p(i), residual_nop(i),residual_p(i), t_nop(i), t_p2(i), t_px(i)] = customGMRES(D, E, b, c, num_iterations(i), precond, a); 
           
            tic;
            Ax = [D.*x_custom(1:length(D),i) + E'*x_custom(length(D)+1:end,i); E*x_custom(1:length(D),i)];
        
            residual_px(i) = norm(Ax-b_tilde)/b_norm; % compute the residual
            
            t_temp = toc;
            t_px(i)= t_px(i)+t_temp;
        else
            [x_custom(:, i), q_nop(i),~, residual_nop(i),~, t_nop(i), ~, ~] = customGMRES(D, E, b, c, num_iterations(i), precond, a); 
        
        end

        if strcmp(mode, 'rate') && i > 1
            rate(i-1)=residual_nop(i)/residual(i-1);
        end
        
        if strcmp(mode, 'minres') 
            fprintf(fileID, "q = %d\n",q_nop(i));
            fprintf(fileID, "r = %d\n",residual_nop(i));
            fprintf(fileID,"t = %d\n", t_nop(i));
            fprintf(fileID,"m = %d\n", num_iterations(i));
        end
 
        if i > 1
            
            del_nop = abs(residual_nop(i)-residual_nop(i-1));
            if del_nop < delta && flag
                display("CustomGMRES stopped since the delta was close to 0 (nop)");
                flag = false;
                position = i;
            end
        end        
    end 
   
    
    
    switch mode
        case 'precond'
            if flag
                position = i;
            end
            scatter(q_nop(1:position),t_nop(1:position),80,"red", "filled");
            hold on;
            scatter(q_p(1:i),t_p2(1:i), 80,"magenta","filled");
            h(1) = scatter(residual_nop(1:position), t_nop(1:position), 35, 'blu','square', 'filled'); % no precond
            h(2) = scatter(residual_p(1:i), t_p2(1:i), 35, 'green', 'filled'); % precond transformed
            h(3) = scatter(residual_px(1:i),t_px(1:i), 45, [0.4940 0.1840 0.5560],'filled'); % precond original
            xlabel('Relative Residual');
            ylabel('Time (s)');
            set(gca, 'XScale', 'log');

            legend(h, ["Non Preconditioned","Preconditioned Transformed","Preconditioned Original"]);
            legend('Location', 'northwest');

        case 'rate'
            ratefunc(num_iterations(1:i-1),rate(1:i-1))
        case 'residual'
            residualfunc(num_iterations(1:i), residual(1:i))
        case 'minres'
            c = 1:i;
            h(1)=scatter(q_nop(1:i),t_nop(1:i),80,"red", "filled");
            hold on;
            h(2)=scatter(residual_nop(1:i), t_nop(1:i), 25, c, 'filled');
            colormap(winter);
            r_minres = zeros(i,1);
            t_minres = zeros(i,1);
            
            fprintf(fileID, "MINRES\n");
            
            
            for j = 1:i
                tic;
                [~,~,r_minres(j)] = minres(A, b_tilde, 1e-15, num_iterations(j));
                t_minres(j) = toc;
                
                fprintf(fileID, "r_minres = %d\n",r_minres(j));
                fprintf(fileID,"t_minsres = %d\n", t_minres(j));
            end
            
            h(3)=scatter(r_minres, t_minres, 60, [0.3 0.2 0.6], 'square', 'filled');
            set(gca, 'XScale', 'log');
            legend(h, ["Lower Bound","CustomGMRES","MINRES"]);
    end

 
