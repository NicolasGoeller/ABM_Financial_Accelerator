function [nw] = network_worth(term, network, actor, recipient)
%Function returns the a vector of the summed value per actor for each
%recipient that is connected to the sctor in the network
% term: a vector (or linear combination of values
% network: a vector of values indicating network connections
% actor: the agents that are connected into
% recipient: the agents that are summed up
nw = zeros(actor, 1);
for i = 1:actor
    I = zeros(1, recipient);
    for  j = 1:recipient
        if i == network(j)
            I(j) = 1;
        else 
            I(j) = 0;
        end
    end
    nw(i) = sum((term .* I));
end
   
    
