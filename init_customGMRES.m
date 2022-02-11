function init_customGMRES(b, c, D, E)


    % given b and c we build b tilde
    b_tilde = [b; c];
    % compute the norm of b tilde
    b_norm = norm(b_tilde);

    % minimum improvements to continue the algorithm
    delta = 10^-8;

    n = length(b_tilde);

    % Now we start to compute the customized GMRES with different number
    % of iterations of the Arnoldi process so that we show how the residual
    % changes with the number of iterations. 

    % different number of iterations of the Arnoldi process proportional to
    % the dimensions of the matrix A
    num_iterations = [50:round(sqrt(n)):round(n/3)];
    

    % initialization of the variables
    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES for different number of iterations 

    t_custom = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES for different number of iterations

    residual = zeros(length(num_iterations), 1); % initialization of the
    % array of the residuals for different number of iterations

    fileid = fopen("test.txt", "w");

    for i=1:length(num_iterations)
        tic; % start counting the time spent for one iteration
        [x_custom(:, i),q] = customGMRES(D, E, b, c, num_iterations(i)); 
        % apply the customized GMRES to the problem
        t_custom(i) = toc; % save the time spent for the customized GMRES 
        % where the number of iterations is num_iterations(i)
        Ax = [D.*x_custom(1:length(D),i)+E'*x_custom(length(D)+1:end,i); 
            E*x_custom(1:length(D),i)];
        residual(i) = norm(Ax-b_tilde)/b_norm; % compute 
        % the residual
        scatter(q, t_custom(i), 80, "red", "filled");
        fprintf(fileid, "m: %d\n", num_iterations(i));
        fprintf(fileid, "q: %d\n", q);
        fprintf(fileid, "residual: %d\n", residual(i));
        fprintf(fileid, "t: %d\n", t_custom(i));
        hold on
        if i>1 && abs(residual(i)-residual(i-1)) < delta
            fprintf("CustomGMRES stopped at the iteration %d because " + ...
                "the difference of the residual is %d\n", i, abs(residual(i)-residual(i-1)));
            break;
        end 
        
    end 
    

    % plotting of the results obtained with customized GMRES
    c = 1:i;
    scatter(residual(1:i), t_custom(1:i), 30, c,'filled');
    colormap(winter);
    hold on


    % We write the labels on the x axis and y axis
    xlabel('Relative residual');
    ylabel('Time (s)');
    
    % We set the values on the x axis 
    % xlim([0,0.3])

    fclose(fileid);
