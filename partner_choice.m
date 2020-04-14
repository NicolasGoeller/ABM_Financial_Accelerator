function [pc] = partner_choice(lender_worth, borrower_leverage, M, network, para)
%Function return a vector of network positions that are selected based on
%the offered prices. There is a tendency to stick with the current partner.
% lender_worth: vector of net_worth of lenders
% borrower_leverage: vector of leverage of borrowers
% M: market visibility share
% network: connection network between lenders and borrowers
% para: lambda in simulationa model
pc = zeros(1, length(borrower_leverage));
for i = 1:length(borrower_leverage)
    sub = zeros(3, round(M*length(lender_worth)));
    sub(1,:) = randsample(length(lender_worth), round(M*length(lender_worth)));
    sub(2,:) = lender_worth(sub(1,:));
    sub(3,:) = sub(2,:) .* borrower_leverage(i);
    
    [cand, I] = min(sub(3,:));
    curr = lender_worth(network(i)) * borrower_leverage(i);
    ps = 0;
    
    if cand < curr
        ps = 1 - exp((para*(cand-curr)/cand));
    end
    
    if rand() < ps
        pc(i) = sub(1,I);
    else
        pc(i) = network(i);
    end
end