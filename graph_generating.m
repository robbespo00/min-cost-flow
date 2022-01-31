function E = graph_generating(n, e)
    
    G=digraph(true(n), 'omitselfloops');
    p=randperm(numedges(G), e);
    G=digraph(G.Edges(p,:));
    E=incidence(G);
    E=E.*round(20*(rand(n,e)-1/2*true(n,e)));
    