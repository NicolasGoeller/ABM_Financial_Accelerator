function [pc] = partner_choice(lender_worth, borrower_leverage, network, borrower_bankrupt, M, lambda, alpha)
%Function return a vector of network positions that are selected based on
%the offered prices. There is a tendency to stick with the current partner.
% lender_worth: vector of net_worth of lenders
% borrower_leverage: vector of leverage of borrowers
% network: connection network between lenders and borrowers
% borrower_bankrupt: vector indicating borrower bankruptcy status
% M: market visibility share
% lambda: partner switching propensity
% alpha: interest rate setting parameter
pc = zeros(1, length(borrower_leverage));
for i = 1:length(borrower_leverage)
    sub = zeros(3, round(M*length(lender_worth)));
    sub(1,:) = randsample(length(lender_worth), round(M*length(lender_worth)));
    sub(2,:) = lender_worth(sub(1,:));
    sub(3,:) = ((sub(2,:).^(alpha*-1)) .* alpha) + ((borrower_leverage(i).^alpha) .* alpha);
    
    [cand, I] = min(sub(3,:));
    curr = ((lender_worth(network(i)).^(alpha*-1)) .* alpha) + ((borrower_leverage(i).^alpha) .* alpha);
    ps = 0;
    
    if cand < curr
        ps = 1 - exp((lambda*(cand-curr)/cand));
    end
    
    if rand() < ps
        pc(i) = sub(1,I);
    else
        pc(i) = network(i);
    end
    %Correction for case the borrower just went bankrupt - Randomize partner
    %selection
    if borrower_bankrupt(i) == 1
        pc(i) = randi([1,length(lender_worth)], 1,1);
    end
    %Correction for case the lender just went bankrupt - Randomize partner
    %selection
    if lender_worth(pc(i)) < 0
        pc(i) = randi([1,length(lender_worth)], 1,1);
    end
end