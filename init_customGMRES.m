function init_customGMRES(A, b, c, D, E)

    % given b and c we build b tilde
    b_tilde = [b; c];
    % compute the norm of b tilde
    b_norm = norm(b_tilde);


    % Now we start to compute the customized GMRES with different number
    % of iterations of the Arnoldi process so that we show how the error
    % changes with the number of iterations. Then, we compute the
    % solution using the MATLAB lsqr() function in order to compare the 
    % results (looking at the time spent by each algorithm and at the error).


    % fixed number of iterations of the Arnoldi process in order to apply
    % customized GMRES
    num_iterations = [50 100 200 500];

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

    
    % At this point we compute the same problem, but instead of applying
    % customized GMRES we use the function implemented in MATLAB in order
    % to check the differences between our solution and the lsqr()
    % solution.
    

    % fixed number of iterations of the lsqr without applying
    % the customized GMRES
    iterations_lsqr = [500, 1000, 2000];

    % initialization of the variables
    x_lsqr = zeros(length(b_tilde), length(iterations_lsqr)); % solutions 
    % of lsqr with different number of iterations
    t_lsqr = zeros(length(iterations_lsqr), 1); % time spent by lsqr with 
    % different number of iterations
    norm_lsqr = zeros(length(iterations_lsqr),1); % norm of (Ax-b) over norm
    % of b, where x is x_lsqr(j) (result given by lsqr function)
    

    for j=1:length(iterations_lsqr)
        tic; % start counting the time spent for one iteration
        x_lsqr(:, j) = lsqr(A, b_tilde, 1e-15, iterations_lsqr(j)); % apply
        % the lsqr to the problem
        t_lsqr(j) = toc; % save the time spent for the lsqr where the 
        % number of iterations is iterations_lsqr(j)
        t_lsqr(j) = round(t_lsqr(j), 2); % ceil the number
        disp(t_lsqr);
        norm_lsqr(j) = norm(A*x_lsqr(:, j)-b_tilde)/b_norm; % compute the norm 
    end

 
    % plotting of the results obtained with lsqr function
    scatter(norm_lsqr, t_lsqr, 50, 'red', 'filled');
    


    % We write the labels on the x axis and y axis
    xlabel('Relative residual');
    ylabel('Time (s)');
    
    % We set the values on the x axis from 0 to 0.15
    xlim([0,0.15])


   
