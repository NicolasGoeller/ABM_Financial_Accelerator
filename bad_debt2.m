function [bd] = bad_debt2(net_worth, credits, network, lender)
% This function returns a vecto of values that represent the summed
% outstanding debt for connected network partners that went bankrupt.
% net_worth:
% credits:
% network: 
% lender
bd = zeros(1,lender);

for i = 1:lender
    I = zeros(1,length(net_worth));
    for  j = 1:length(net_worth)
        if i == network(j) && net_worth(j) <= 0
            I(j) = 1;
        else 
            I(j) = 0;
        end
    end
    
    bd(i) = sum(credits .* I);
end