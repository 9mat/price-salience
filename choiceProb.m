function [ prob ] = choiceProb( theta, Data, n )
%CHOICEPROB Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta, n);

if n.demogr > 0
    alphai = params.alpha + Data.demogr*params.gamma;
else
    alphai = params.alpha;
end
V = bsxfun(@times, alphai, Data.price) + Data.X*params.beta;

for choice = 2:n.choice
    index = Data.choice == choice;
    V(index, :) = V(index, :)*n.M(:,:,choice)';
end

prob = zeros(size(V,1),1);

for choice = 1:n.choice
    for treat = 1:n.treat
        index = (Data.choice == choice) & (Data.treat == treat);
        S = n.M(:,:,choice)*params.S(:,:,treat);
        try
            prob(index) = mvncdf(-V(index,:), zeros(1, n.choice-1), S*S');
        catch
            fprintf('error: probably singular S*S\n');
            prob(index) = 0;
        end
    end
end

end

