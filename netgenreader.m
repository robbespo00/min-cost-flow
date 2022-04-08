function [E, b, c] = netgenreader(filename)

    
    M = readmatrix(filename);
    
    c = M(1:1024,3);
    
    G = digraph(M(1025:end,2), M(1025:end,3), 'omitselfloops');
    E = incidence(G);
    
    
    b = M(1025:end, 4);