function init_customGMRES(filename)
    
    M = readmatrix(filename);
    G = digraph(M(:,1), M(:,2), 'omitselfloops');
    E = incidence(G);
    E = M(:,3)'.*E;
    [nodi, e] = size(E);
    
    [D, b, c] = matrices_generator(nodi, e, 14, 17, 14);
    A = [D.*eye(e) E'; E zeros(nodi)];

    % given b and c we build b tilde
    b_tilde = [b; c];
    % compute the norm of b tilde
    b_norm = norm(b_tilde);
    n = nodi + e;
    delta = 10^(-10);


    % Now we start to compute the customized GMRES with different number
    % of iterations of the Arnoldi process so that we show how the error
    % changes with the number of iterations. Then, we compute the
    % solution using the MATLAB lsqr() function in order to compare the 
    % results (looking at the time spent by each algorithm and at the error).


    % fixed number of iterations of the Arnoldi process in order to apply
    % customized GMRES
    num_iterations = [50:round(sqrt(n)):round(n/15)];

    % initialization of the variables
    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES with different number of iterations 
    t_custom = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES with different number of iterations
    residual = zeros(length(num_iterations), 1); % norm of (Ax-b) over 
    % norm of b, where x is the result given by customized GMRES

    for i = 1:length(num_iterations)
        tic; % start counting the time spent for one iteration
        [x_custom(:, i),q] = customGMRES(D, E, b, c, num_iterations(i)); 
        % apply the customized GMRES to the problem
        t_custom(i) = toc; % save the time spent for the customized GMRES 
        % where the number of iterations is num_iterations(j)
        t_custom(i) = round(t_custom(i), 2); % approximate the time up to 2 
        % digits after the comma
        Ax = [D.*x_custom(1:length(D),i) + E'*x_custom(length(D)+1:end,i); E*x_custom(1:length(D),i)];
        residual(i) = norm(Ax-b_tilde)/b_norm; % compute 
        % the residual

        scatter(q,t_custom(i),70,"red", "filled");
        hold on
        if i > 1
            del = abs(residual(i)-residual(i-1));
            if del < delta
                break;
            end
        end

        fprintf("q = %d\n",q);
        fprintf("r = %d\n",residual(i));
        fprintf("t = %d\n", t_custom(i));
        fprintf("m = %d\n", num_iterations(i));

    end 


    % plotting of the results obtained with customized GMRES
    c = 1:i;
    scatter(residual(1:i), t_custom(1:i), 25, c,'filled');
    colormap(winter);
    hold on


    % We write the labels on the x axis and y axis
    xlabel('Relative residual');
    ylabel('Time (s)');
    
    % We set the values on the x axis from 0 to 0.15
    xlim([0,0.3]);

    fprintf("delta = %d", del);

    r_minres = zeros(i,1);
    t_minres = zeros(i,1);

    for j = 1:i
        tic;
        [~,~,r_minres(j)] = minres(A, b_tilde, 1e-6, num_iterations(j));
        t_minres(j) = toc;

        fprintf("r = %d\n",r_minres(j));
        fprintf("t = %d\n", t_minres(j));
    end

    scatter(r_minres, t_minres, 80, c, 'p', 'filled');

    
   
    
