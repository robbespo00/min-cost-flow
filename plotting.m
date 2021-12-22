function plotting(A, x, x2, b, t1, t2)
    % We save the norm of b to compute the relative residual
    b_norm = norm(b);
    % We compute the norm of the solution given by our algorithm
    normQR = norm(A*x-b)/b_norm;
    % We compute the norm of the solution given by the standard least
    % square
    standard_norm = norm(A*x2-b)/b_norm;
    % We plot our solution where the x axis is the relative residual and
    % the y axis is the computational time
    scatter(normQR, t1, 50, 'filled', 'green');
    hold on;
    % We plot standard solution where the x axis is the relative residual and
    % the y axis is the computational time
    scatter(standard_norm, t2, 50, 'filled', 'red');
    % We set the x axis from 0 to 0.1
    xlim([0,0.1])
    % We print the legend of the plot
    legend('Our algorithm', 'Standard LSQR')
    % We ceil the value of the norm
    our_norm_standard = sprintf("%g", round(normQR,4));
    % We ceil the value of the norm
    standard_norm_standard = sprintf("%g", round(standard_norm,4));
    % We write the coordinates of the points
    text_1 ="   ("+ our_norm_standard + ", " + t1 + ")";
    text_2 ="   ("+ standard_norm_standard + ", " + t2 + ")";
    text(normQR, t1,text_1);
    text(standard_norm, t2,text_2);
    % We write the labels on the x axis and y axis
    xlabel('Relative residual');
    ylabel('Time (s)');
    