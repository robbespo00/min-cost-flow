function [E, b, c] = netgenreader(filename)

    
    M = readmatrix(filename);
    
    c = M(1:4096,3);
    
    G = digraph(M(4097:end,2), M(4097:end,3), 'omitselfloops');
    E = incidence(G);
    b = M(4097:end, 4);
    