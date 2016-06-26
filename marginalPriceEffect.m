function [ ame ] = marginalPriceEffect( theta, Data, n )
%MARGINALPRICEEFFECT Summary of this function goes here
%   ame(:, i, j): marginal effect of changing price of i on j market share
prob0 = zeros(n.obs, n.choice, n.choice);
prob1 = zeros(n.obs, n.choice, n.choice);
d = 1e-5;

for i = 1:n.choice
    if i == 1 
        Data.price = Data.price + d;
    else
        Data.price(:, i-1) = Data.price(:, i-1) - d;
    end
    
    for j = 1:n.choice
        Data.choice(:) = j;
        prob0(:,j,i) = choiceProb(theta, Data, n);
    end
    
    if i == 1 
        Data.price = Data.price - 2*d;
    else
        Data.price(:, i-1) = Data.price(:, i-1) + 2*d;
    end
    
    for j = 1:n.choice
        Data.choice(:) = j;
        prob1(:,j,i) = choiceProb(theta, Data, n);
    end
    
    if i == 1 
        Data.price = Data.price + d;
    else
        Data.price(:, i-1) = Data.price(:, i-1) - d;
    end

end

ame = (prob1 - prob0)/(2*d);

ame = squeeze(mean(ame, 1));

end

