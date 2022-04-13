function [E, b, c] = netgenreader(filename, nodes)

    
    M = readmatrix(filename);
    
    c = M(1:nodes,3);
    
    G = digraph(M(nodes+1:end,2), M(nodes+1:end,3), 'omitselfloops');
    E = incidence(G);
    
    
    b = M(nodes+1:end, 4);