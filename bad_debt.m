function [bd] = bad_debt(net_worth, credits, network, lender)
% This function returns a vector of values that represent the summed
% outstanding debt for connected network partners that went bankrupt.
% net_worth:
% credits:
% network: 
% lender
bd = zeros(1,lender);
for i = 1:lender
    bd(i) = sum(credits((net_worth <= 0) & (network == i)));
end
