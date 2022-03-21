function init_customGMRES(filename, mode, generate, distribution)
    
    [E, b, c] = netgenreader(filename);
    [nodi, e] = size(E);
    
    % given b and c we build b tilde
    b=-b;
    b_tilde = [b; c];
    
    % compute the norm of b tilde
    b_norm = norm(b_tilde);
    n = nodi + e;
    delta = 10^(-10);
    
    
    
    
    
    if generate
        rng(461);
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
            jump = sqrt(n);
        case 'rate'
            jump=1;
        case 'residual'
            jump=1;
    end

    
    
    num_iterations = [50:jump:round(n/32)];

    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES with different number of iterations 
    t_custom = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES with different number of iterations
    residual = zeros(length(num_iterations), 1); 
    rate = zeros(length(num_iterations)-1, 1);
    
    for i = 1:length(num_iterations)
        disp(i);
        tic; % start counting the time spent for one iteration
        [x_custom(:, i), q] = customGMRES(D, E, b, c, num_iterations(i)); 
        % apply the customized GMRES to the problem
        t_custom(i) = toc; % save the time spent for the customized GMRES 
        % where the number of iterations is num_iterations(j)
        %t_custom(i) = round(t_custom(i), 2); % approximate the time up to 2 
        % digits after the comma
        Ax = [D.*x_custom(1:length(D),i) + E'*x_custom(length(D)+1:end,i); E*x_custom(1:length(D),i)];
        residual(i) = norm(Ax-b_tilde)/b_norm; % compute 
        % the residual
        
        if strcmp(mode, 'rate') && i > 1
            rate(i-1)=residual(i)/residual(i-1);
        end
        
        if strcmp(mode, 'minres')
            scatter(q,t_custom(i),70,"red", "filled");
            fprintf(fileID, "q = %d\n",q);
            fprintf(fileID, "r = %d\n",residual(i));
            fprintf(fileID,"t = %d\n", t_custom(i));
            fprintf(fileID,"m = %d\n", num_iterations(i));
        end
 
        if i > 1
            del = abs(residual(i)-residual(i-1));
            if del < delta
                display("CustomGMRES stopped since the delta was close to 0");
                break;
            end
        end

        
        
    end 
   
    
    switch mode
        case 'rate'
            scatter(num_iterations(1:i-1),rate(1:i-1),'green','filled');
            ratefunc()
        case 'residual'
            scatter(num_iterations(1:i), residual(1:i),'blue', 'filled');
            residualfunc()
        case 'minres'
            c = 1:i;
            scatter(residual(1:i), t_custom(1:i), 25, c, 'filled');
            colormap(winter);
            r_minres = zeros(i,1);
            t_minres = zeros(i,1);
            
            fprintf(fileID, "MINRES\n");
            
            for j = 1:i
                tic;
                [~,~,r_minres(j)] = minres(A, b_tilde, 1e-15, num_iterations(j));
                t_minres(j) = toc;
                
                fprintf(fileID, "r = %d\n",r_minres(j));
                fprintf(fileID,"t = %d\n", t_minres(j));
            end
            
            scatter(t_minres, r_minres, 60, [0.3 0.2 0.6], 'square', 'filled');
            set(gca, 'YScale', 'log');
            set(gca, 'XScale', 'log');
            h(1) = scatter(-1,-1,40,"blue", 'filled');
            h(2) = scatter(-1,-1,80,"blue", 'square','filled');
            legend(h, ["CustomGMRES","MINRES"]);
    end

 
%     dimension = e;
%     t = t_custom(i-1);
%     scatter(e+nodi, t_custom(i-1), 25, 'green', 'filled');
%     scatter(e, residual(1:i), 25, c, 'filled');

    
