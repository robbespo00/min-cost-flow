function [E, D, b, c] = matrices_generator(n, e, d_range, e_range, b_range, c_range)

    E=graph_generator(n,e, e_range);
    D = d_range*rand(e,1);
    b= 2*b_range*(rand(e,1)-(ones(e,1))/2);
    c= 2*c_range*(rand(n,1)-(ones(n,1))/2);
    