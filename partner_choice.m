function [pc] = partner_choice(price, M, recipient, network, para)
%Function return a vector of network positions that are selected based on
%the offered prices. There is a tendency to stick with the current partner.
% price: vector of prices offered up in the round
% recipient
pc = zeros(1, recipient);
for i = 1:recipient
    sub = zeros(2, M*length(price));
    sub(1,:) = randsample(length(price), M*length(price));
    sub(2,:) = price(sub(1,:));
    
    [cand, I] = min(sub(2,:));
    ps = 0;
    if cand < price(i)
        ps = 1 - exp((para*(cand-price(i))/cand));
    end
    
    if rand() < ps
        pc(i) = sub(1,I);
    else
        pc(i) = network(i);
    end
end