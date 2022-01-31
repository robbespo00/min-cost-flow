function init_customGMRES(A, b, c, D, E)

    % given b and c we build b tilde
    b_tilde = [b; c];
    % compute the norm of b tilde
    b_norm = norm(b_tilde);

    n = length(A);


    % Now we start to compute the customized GMRES with different number
    % of iterations of the Arnoldi process so that we show how the error
    % changes with the number of iterations. Then, we compute the
    % solution using the MATLAB lsqr() function in order to compare the 
    % results (looking at the time spent by each algorithm and at the error).


    % fixed number of iterations of the Arnoldi process in order to apply
    % customized GMRES
    num_iterations = [50:round(sqrt(n)):round(n/4), round(n/2), n];

    % initialization of the variables
    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES with different number of iterations 
    t_custom = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES with different number of iterations
    norm_custom = zeros(length(num_iterations), 1); % norm of (Ax-b) over 
    % norm of b, where x is the result given by customized GMRES

    for i=1:length(num_iterations)
        tic; % start counting the time spent for one iteration
        x_custom(:, i) = customGMRES(D, E, b, c, num_iterations(i)); 
        % apply the customized GMRES to the problem
        t_custom(i) = toc; % save the time spent for the customized GMRES 
        % where the number of iterations is num_iterations(j)
        t_custom(i) = round(t_custom(i), 2); % approximate the time up to 2 
        % digits after the comma
        norm_custom(i) = norm(A*x_custom(:,i)-b_tilde)/b_norm; % compute 
        % the norm
    end 


    % plotting of the results obtained with customized GMRES
    c = 1:length(num_iterations);
    scatter(norm_custom, t_custom, 50, c,'filled');
    colormap(winter);
    hold on


    % We write the labels on the x axis and y axis
    xlabel('Relative residual');
    ylabel('Time (s)');
    
    % We set the values on the x axis from 0 to 0.15
    xlim([0,0.3])
