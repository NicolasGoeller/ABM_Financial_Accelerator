function [nw] = network_worth(term, network, actor)
%Function returns the a vector of the summed value per actor for each
%recipient that is connected to the sctor in the network
% term: a vector (or linear combination of values
% network: a vector of values indicating network connections
% actor: the agents that are connected into
nw = zeros(1, actor);
for i = 1:actor
    nw(i) = sum(term(network == i));
end
   
    
