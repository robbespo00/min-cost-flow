function compl_dim(filename, num_iterations)


  
   dimension = zeros(num_iterations+1, 1);
   t = zeros(num_iterations+1, 1);
 
    for i=0:num_iterations
       
       i
      
       net = strcat(strcat(erase(filename, ".txt"),int2str(i)),'.txt');
       
       [dimension(i+1), t(i+1)] = init_customGMRES(net);
       
    end
    scatter(dimension, t, 35,[0.9290 0.6940 0.1250], 'filled');
    %[244/255, 110/255, 6/255]
    p = polyfit(dimension,t,1);
    fplot(@(t) p(1)*t+p(2), 'red', 'LineWidth', 1.7);