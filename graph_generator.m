function E = graph_generator(n, e, e_range)
    
    G=digraph(true(n), 'omitselfloops');
    p=randperm(numedges(G), e);
    G=digraph(G.Edges(p,:));
    E=incidence(G);
    E=E.*round(2*e_range*(rand(n,e)-1/2*true(n,e)));
    