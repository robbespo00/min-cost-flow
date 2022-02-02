function init_customGMRES(A, b, c, D, E)

    % We give as input the matrix A in order to compute the residual, this
    % does not influence the real computation time since the calculation is 
    % outside the "tic toc" function
    n = length(A);

    % given b and c we build b tilde
    b_tilde = [b; c];
    % compute the norm of b tilde
    b_norm = norm(b_tilde);


    % Now we start to compute the customized GMRES with different number
    % of iterations of the Arnoldi process so that we show how the residual
    % changes with the number of iterations. 

    % different number of iterations of the Arnoldi process proportional to
    % the dimensions of the matrix A
    num_iterations = [50:round(sqrt(n)):n-1, n];
    

    % initialization of the variables
    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES for different number of iterations 

    t_custom = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES for different number of iterations

    norm_custom = zeros(length(num_iterations), 1); % initialization of the
    % array of the residuals for different number of iterations


    for i=1:length(num_iterations)
        tic; % start counting the time spent for one iteration
        x_custom(:, i) = customGMRES(D, E, b, c, num_iterations(i)); 
        % apply the customized GMRES to the problem
        t_custom(i) = toc; % save the time spent for the customized GMRES 
        % where the number of iterations is num_iterations(i)
        norm_custom(i) = norm(A*x_custom(:,i)-b_tilde)/b_norm; % compute 
        % the residual
    end 


    % plotting of the results obtained with customized GMRES
    c = 1:length(num_iterations);
    scatter(norm_custom, t_custom, 50, c,'filled');
    colormap(winter);
    hold on


    % We write the labels on the x axis and y axis
    xlabel('Relative residual');
    ylabel('Time (s)');
    
    % We set the values on the x axis 
    xlim([0,0.3])
