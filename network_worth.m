function nw = network_worth(term, network, actor, recipient)
nw = zeros(1:actor)
for i = 1:actor
    I = zeros(1:recipient)
    for  j = 1:recipient
        if i = network(j)
            I(j) = 1
        else 
            I(j) = 0
    nw(i) = sum(term .* I)
