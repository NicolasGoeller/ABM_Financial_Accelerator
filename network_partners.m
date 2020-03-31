function [np] = network_partners(term, network, recipient)
% Function returns vector of length recipient that contains the value of
% term for the network partner of each column actor in network
% term: a vector 
% network: a vector of values indicating network connections
% recipient: the agents that are connected into
np = zeros(1,recipient);
for i = 1:recipient
    p = network(i);
    np(i) = term(p);
end