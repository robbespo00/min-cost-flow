function [E, b, c] = netgenreader(filename)

    
    M = readmatrix(filename);
    
    c = M(1:256,3);
    
    G = digraph(M(257:end,2), M(257:end,3), 'omitselfloops');
    E = incidence(G);
    
    
    b = M(257:end, 4);
    